#include <filesystem>
#include <fstream>
#include <iomanip>
#include <set>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/mesh/euler/eu_prim.h>

#include <zephyr/geom/indexing.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/geom/grid.h>

#include <zephyr/mesh/amr/apply.h>
#include <zephyr/mesh/amr/balancing.h>
#include <zephyr/mesh/amr/rotations.h>

#include <zephyr/utils/json.h>
#include <zephyr/io/pvd_file.h>

namespace zephyr::mesh {

using namespace geom;
using namespace utils;
using generator::Rectangle;
using generator::Cuboid;
using utils::Stopwatch;
using namespace io;
using namespace geom::indexing;

namespace {

// Генерация локальной сетки из grid
AmrCells make_locals(geom::Grid&& grid, bool unique_nodes) {
    // Опции генерации сетки
    constexpr Grid::BuildOptions options{
        .build_faces=true
    };

    grid.finalize(options);

    if (!grid.has_faces()) {
        throw std::runtime_error("Grid was built with wrong options, need faces per cell");
    }

    AmrCells locals({
        .dim = grid.dimension(),
        .adaptive = grid.adaptive(),
        .linear = true,
        .axial = false,
        .nodes = unique_nodes
    });

    if (!locals.adaptive()) {
        locals.resize(grid.n_cells(), grid.total_faces_per_cell(), grid.total_nodes_per_cell());
    }
    else {
        locals.resize_amr(grid.n_cells());
    }

    const auto& nodes = grid.nodes();
    const auto& cells = grid.cells();

    locals.faces.offsets[0] = 0;
    locals.verts.offsets[0] = 0;
    for (index_t ic = 0; ic < locals.n_cells(); ++ic) {
        locals.rank[ic] = 0;
        locals.index[ic] = ic;

        locals.flag[ic] = 0;
        locals.b_idx[ic] = ic;
        locals.z_idx[ic] = 0;
        locals.level[ic] = 0;

        locals.volume[ic] = cells[ic].volume();
        locals.center[ic] = cells[ic].centroid();

        if (locals.axial()) {
            // locals.volume_alt[ic] = geom.cell_volumes_alt[ic];
        }

        // Число узлов и граней ячейки
        int n_nodes = cells[ic].n_nodes();
        int n_faces = cells[ic].n_faces();
        int max_faces = n_faces;
        if (locals.adaptive()) {
            max_faces = locals.dim() < 3 ? Side2D::n_subfaces() : Side3D::n_subfaces();
        }

        // Выставить индексы грани
        locals.faces.offsets[ic + 1] = locals.faces.offsets[ic] + max_faces;
        if (cells[ic].type() != CellType::POLYHEDRON) {
            locals.faces.insert(locals.faces.offsets[ic], cells[ic].type(), max_faces);
        }
        else {
            auto iface = locals.faces.offsets[ic];
            for (int i = 0; i < n_faces; ++i) {
                const auto& face = cells[ic].get_face(i);
                locals.faces.vertices[iface + i].fill(-1);
                for (int j = 0; j < face.n_nodes(); ++j) {
                    locals.faces.vertices[iface + i][j] = face.node_idx(j);
                }
                locals.faces.set_undefined(iface + i);
            }
        }
        const auto& node_ids = cells[ic].nodes();

        locals.verts.offsets[ic + 1] = locals.verts.offsets[ic] + n_nodes;
        for (int i = 0; i < n_nodes; ++i) {
            locals.verts[locals.verts.offsets[ic] + i] = nodes[node_ids[i]].pos;
        }

        // Геометрия
        for (int i = 0; i < n_faces; ++i) {
            const auto& face = cells[ic].get_face(i);

            int iface = locals.faces.offsets[ic] + i;
            locals.faces.boundary[iface] = face.bc();
            locals.faces.area[iface]     = face.area();
            locals.faces.center[iface]   = face.center();
            locals.faces.normal[iface]   = face.normal();

            if (locals.axial()) {
                //locals.faces.area_alt[iface] = grid_geom.face_areas_alt[jface];
            }
        }

        // Смежность
        for (int i = 0; i < n_faces; ++i) {
            const auto& face = cells[ic].get_face(i);

            auto iface = locals.faces.offsets[ic] + i;

            locals.faces.adjacent.rank[iface]  = 0;
            locals.faces.adjacent.index[iface] = face.neib();
            locals.faces.adjacent.ghost[iface] = -1;
            locals.faces.adjacent.basic[iface] = ic;
            locals.faces.adjacent.rotation[iface] = 0;

            if (locals.faces.is_boundary(iface)) {
                locals.faces.adjacent.index[iface] = ic;
            }
        }
    }
    return locals;
}

} // anonymous namespace

void EuMesh::sync_params_() {
#ifdef ZEPHYR_MPI
    if (!mpi::single()) {
        int dim = local_cells_.dim();
        int adapt = local_cells_.adaptive();
        int axial = local_cells_.axial();
        int nodes = local_cells_.verts.has_nodes();

        // Соберем базовые характеристики сетки с master-процесса
        mpi::broadcast(0, dim);
        mpi::broadcast(0, adapt);
        mpi::broadcast(0, axial);
        mpi::broadcast(0, nodes);

        if (!mpi::master()) {
            MeshOpts opts {
                .dim      = dim,
                .adaptive = bool(adapt),
                .axial    = bool(axial),
                .nodes    = bool(nodes)
            };
            local_cells_ = AmrCells(opts);
        }
    }

    // установить опции
    tourists_.init_types(local_cells_);
    migrants_.init_types(local_cells_);
#endif
}

void EuMesh::build_(Generator& gen, bool unique_nodes) {
    if (mpi::master()) {
        if (gen.can_make_cells()) {
            local_cells_ = gen.make_cells(unique_nodes);
        }
        else {
            Grid grid = gen.make();
            local_cells_ = make_locals(std::move(grid), unique_nodes);
        }
        init_amr();
        if (unique_nodes) {
            local_nodes_.setup_for(local_cells_);
        }
    }
    sync_params_();

    // Иногда генератор позволяет восстановить структуру сетки
    if (mpi::single()) {
        // Структурированность не работает для распределенных расчетов
        if (gen.can_cast<Rectangle>()) {
            Rectangle& rect = gen.cast<Rectangle>();
            if (rect.structured()) {
                structured_ = true;
                nx_ = rect.nx();
                ny_ = rect.ny();
                nz_ = 1;
            }
        } else if (gen.can_cast<Cuboid>()) {
            Cuboid& cuboid = gen.cast<Cuboid>();
            structured_ = true;
            nx_ = cuboid.nx();
            ny_ = cuboid.ny();
            nz_ = cuboid.nz();
        }
    }
}

EuMesh::EuMesh(Grid&& grid, bool unique_nodes) {
    if (mpi::master()) {
        local_cells_ = make_locals(std::move(grid), unique_nodes);
        init_amr();
        if (unique_nodes) {
            local_nodes_.setup_for(local_cells_);
        }
    }
    sync_params_();
}

EuMesh::EuMesh(Generator& gen, bool unique_nodes) {
    build_(gen, unique_nodes);
}

EuMesh::EuMesh(const Json& config) {
    // Создать сетку на master-процессе
    auto gen = Generator::create(config);

    bool adaptive = false;
    if (config["adaptive"]) {
        adaptive = config["adaptive"].as<bool>();
    }
    bool unique_nodes = false;
    if (config["nodes"]) {
        unique_nodes = config["nodes"].as<bool>();
    }
    max_level_ = 0;
    if (config["max_level"]) {
        max_level_ = std::max(0, config["max_level"].as<int>());
        if (max_level_ > 15) {
            std::cerr << "Max level set up to " << max_level_ << ", decreased to 15\n";
            max_level_ = 15;
        }
        // Есть ключевое слово adaptive, выставлено на false
        if (!adaptive) {
            max_level_ = 0;
        }
    }

    build_(*gen, unique_nodes);

    // Если есть декомпозиция, то использовать
    if (!mpi::single()) {
        Box domain = bbox();

        if (config["decomp"]) {
            // Там все свойства декомпозиции автоматически ставятся
            decomp_ = Decomposition::create(domain, config["decomp"]);

            //set_decomposition("XY");
            set_decomposition(decomp_, true);
        }
        else {
            // Определяем размерность сетки
            int dim = local_cells_.empty() ? 0 : local_cells_.dim();
            dim = mpi::max(dim);

            z_assert((dim == 2 || dim == 3), "Strange dimension, EuMesh constructed by json");

            // По умолчанию что-то такое
            set_decomposition(dim < 3 ? "XY" : "XYZ");
        }
    }
}

EuMesh EuMesh::PolySet(int dim) {
    EuMesh mesh;
    MeshOpts opts{.dim=dim, .adaptive=false, .nodes=false};
    mesh.local_cells_ = AmrCells(opts);
    return mesh;
}

void EuMesh::init_amr() {
    if (!local_cells_.adaptive()) return;

    amr::find_rotations(local_cells_);
}

bool EuMesh::adaptive() const {
    return max_level_ > 0;
}

int EuMesh::max_level() const {
    return max_level_;
}

void EuMesh::set_max_level(int max_level) {
    if (!local_cells_.adaptive()) {
        max_level_ = 0;
    } else {
        max_level_ = std::max(0, std::min(max_level, 15));
    }
}

void EuMesh::set_distributor(const std::string& name) {
    if (name == "empty") {
        distributor_ = Distributor::empty();
    } else {
        distributor_ = Distributor::simple();
    }
}

void EuMesh::set_distributor(Distributor distr) {
    distributor_ = std::move(distr);
}

void EuMesh::balance_flags() {
#if SCRUTINY
    static bool first_time = true;
    if (first_time) {
        if (check_base() < 0) {
            throw std::runtime_error("Check base failed");
        }
        first_time = false;
    }
#endif
    if (mpi::single()) {
        amr::balance_flags(local_cells_, max_level_);
    }
#ifdef ZEPHYR_MPI
    else {
        amr::balance_flags(local_cells_, max_level_, tourists_);
    }
#endif
}

void EuMesh::apply_flags() {
#if SCRUTINY
    static size_t pvd_counter = 0;

    static PvdFile locals_before("app_bef_locals", "debug");
    static PvdFile ghosts_before("app_bef_ghosts", "debug");
    static PvdFile border_before("app_bef_border", "debug");
    static PvdFile locals_after("app_aft_locals", "debug");
    static PvdFile ghosts_after("app_aft_ghosts", "debug");
    static PvdFile border_after("app_aft_border", "debug");

    if (pvd_counter == 0) {
        Variables vars = {"rank", "index", "next", "level", "flag", "faces2D"};
        locals_before.variables = vars;
        ghosts_before.variables = vars;
        border_before.variables = vars;
        locals_after.variables = vars;
        ghosts_after.variables = vars;
        border_after.variables = vars;
    }
    //locals_before.save(local_cells_, pvd_counter);
    //ghosts_before.save(ghost_cells_, pvd_counter);
    //border_before.save(tourists_.border_cells_, pvd_counter);
    mpi::barrier();
#endif

    if (mpi::single()) {
        amr::apply(local_cells_, distributor_);
    }
#ifdef ZEPHYR_MPI
    else {
        amr::apply(local_cells_, distributor_, tourists_);
    }
#endif

#if SCRUTINY
    //locals_after.save(local_cells_, pvd_counter);
    //ghosts_after.save(ghost_cells_, pvd_counter);
    //border_after.save(tourists_.border_cells_, pvd_counter);
    mpi::barrier();
    ++pvd_counter;

    mpi::for_each([&]() {
        if (check_refined() < 0) {
            std::cout << "Check refined for rank " << mpi::rank() << "\n";
            throw std::runtime_error("Check refined failed");
        }
    });
#endif
}

void EuMesh::make_shuba(int count) {
    count = std::max(0, std::min(count, 50));

    for (int i = 1; i <= count; i++) {
#ifdef ZEPHYR_MPI
        if (!mpi::single()) {
            tourists_.sync<MpiTag::FLAG>(local_cells_);
        }
#endif
        for_each([&](EuCell& cell) {
            if (cell.flag() > 0) return;
            for (auto face: cell.faces()) {
                if (face.neib_flag() == i) {
                    cell.set_flag(i + 1);
                    return;
                }
            }
        });
    }
}

void EuMesh::refine() {
    if (!adaptive()) { return; }
    structured_ = false;

    static Stopwatch balance;
    static Stopwatch apply;
    static Stopwatch full;

    // Для однопроцессорной версии при пустой сетке сразу выход
    if (mpi::single() && local_cells_.empty()) {
        throw std::runtime_error("EuMesh::refine(): Empty mesh");
    }

    if (local_cells_.has_nodes()) {
        throw std::runtime_error("EuMesh::refine(): Unique nodes are not supported");
    }

    full.resume();

    balance.resume();
    balance_flags();
    balance.stop();

    apply.resume();
    apply_flags();
    apply.stop();

    full.stop();


    /*//R version
    threads::parallel_for(0, local_cells_.size(),
        [&locals=local_cells_](index_t ic) {
            index_t iface = locals.faces.offsets[ic];
            std::array<double, 4> areas;
            for (Side2D side: Side2D::items()) {
                areas[side] = locals.faces.area[iface + side];
                if (locals.faces.is_actual(iface + side[1])) {
                    areas[side] += locals.faces.area[iface + side[1]];
                }
            }
            double avg_lr = 0.5 * (areas[Side2D::L] + areas[Side2D::R]);
            double avg_bt = 0.5 * (areas[Side2D::B] + areas[Side2D::T]);
            areas[Side2D::L] = avg_lr;
            areas[Side2D::R] = avg_lr;
            areas[Side2D::B] = avg_bt;
            areas[Side2D::T] = avg_bt;

            for (Side2D side: Side2D::items()) {
                if (locals.faces.is_undefined(iface + side[1])) {
                    locals.faces.area[iface + side] = areas[side];
                }
                else {
                    locals.faces.area[iface + side[0]] = 0.5 * areas[side];
                    locals.faces.area[iface + side[1]] = 0.5 * areas[side];
                }
            }
        });*/
    //Rversion



#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        mpi::cout << "  Balance flags: " << std::setw(14) << balance.milliseconds_mpi() << " ms\n";
        mpi::cout << "  Apply flags:   " << std::setw(14) << apply.milliseconds_mpi() << " ms\n";
        mpi::cout << "Refine elapsed:  " << std::setw(14) << full.milliseconds_mpi() << " ms\n";
    }
    ++counter;
#endif
}

void EuMesh::refine_full(int level) {
    if (level < 0 || level > max_level_) {
        level = max_level_;
    }
    if (level == 0) {
        return;
    }

    for (int i = 0; i < level; ++i) {
        for_each([level](EuCell &cell) {
            cell.set_flag(cell.level() < level ? 1 : 0);
        });
        refine();
    }
}

void EuMesh::check_reference(bool fix) {
    Box box = bbox();
    double xmin = box.vmin.x();
    double xmax = box.vmax.x();
    double ymin = box.vmin.y();
    double ymax = box.vmax.y();

    if (local_cells_.empty()) return;

    int level = local_cells_.level[0];

    int nx = nx_ * std::pow(2, level);
    int ny = ny_ * std::pow(2, level);

    double hx = (xmax - xmin) / nx;
    double hy = (ymax - ymin) / ny;
    double volume = hx * hy;


    std::cout << std::setprecision(14);

    std::cout << "BBox: " << box << "\n";

    std::cout << "Structured: " << nx << " x " << ny << " cells\n";

    auto get_center = [=](int i, int j) -> Vector3d {
        return {xmin + (i + 0.5) * hx, ymin + (j + 0.5) * hy, 0.0};
    };
    auto lface_center = [=](int i, int j) -> Vector3d {
        return {xmin + (i + 0.0) * hx, ymin + (j + 0.5) * hy, 0.0};
    };
    auto rface_center = [=](int i, int j) -> Vector3d {
        return {xmin + (i + 1.0) * hx, ymin + (j + 0.5) * hy, 0.0};
    };
    auto bface_center = [=](int i, int j) -> Vector3d {
        return {xmin + (i + 0.5) * hx, ymin + (j + 0.0) * hy, 0.0};
    };
    auto tface_center = [=](int i, int j) -> Vector3d {
        return {xmin + (i + 0.5) * hx, ymin + (j + 1.0) * hy, 0.0};
    };
    auto get_vertex = [=](double i, double j) -> Vector3d {
        return {xmin + i * hx, ymin + j * hy, 0.0};
    };

    double cell_volume_error = 0.0;
    double cell_center_error = 0.0;
    double face_area_error = 0.0;
    double face_normal_error = 0.0;
    double face_center_error = 0.0;
    double vertices_error = 0.0;

    for (auto& cell: local_cells_) {
        int i = std::round((cell.x() - 0.5 * hx - xmin) / hx);
        int j = std::round((cell.y() - 0.5 * hy - ymin) / hy);

        // Погрешности данных ячейки
        cell_volume_error = std::max(cell_volume_error, std::abs(cell.volume() - volume));
        cell_center_error = std::max(cell_center_error, (cell.center() - get_center(i, j)).norm());

        // Погрешности для граней
        face_area_error = std::max(
                {
                        face_area_error,
                        std::abs(cell.face(Side2D::L).area() - hy),
                        std::abs(cell.face(Side2D::R).area() - hy),
                        std::abs(cell.face(Side2D::B).area() - hx),
                        std::abs(cell.face(Side2D::T).area() - hx),
                });

        face_normal_error = std::max(
                {
                        face_normal_error,
                        (cell.face(Side2D::L).normal() + Vector3d::UnitX()).norm(),
                        (cell.face(Side2D::R).normal() - Vector3d::UnitX()).norm(),
                        (cell.face(Side2D::B).normal() + Vector3d::UnitY()).norm(),
                        (cell.face(Side2D::T).normal() - Vector3d::UnitY()).norm()
                });

        face_center_error = std::max(
                {
                        face_center_error,
                        (cell.face(Side2D::L).center() - lface_center(i, j)).norm(),
                        (cell.face(Side2D::R).center() - rface_center(i, j)).norm(),
                        (cell.face(Side2D::B).center() - bface_center(i, j)).norm(),
                        (cell.face(Side2D::T).center() - tface_center(i, j)).norm()
                });

        // Погрешности для вершин
        SqQuad quad = cell.mapping<2>();
        vertices_error = std::max(
                {
                        vertices_error,
                        (quad.vs<-1, -1>() - get_vertex(i + 0.0, j + 0.0)).norm(),
                        (quad.vs< 0, -1>() - get_vertex(i + 0.5, j + 0.0)).norm(),
                        (quad.vs<+1, -1>() - get_vertex(i + 1.0, j + 0.0)).norm(),
                        (quad.vs<-1,  0>() - get_vertex(i + 0.0, j + 0.5)).norm(),
                        (quad.vs< 0,  0>() - get_vertex(i + 0.5, j + 0.5)).norm(),
                        (quad.vs<+1,  0>() - get_vertex(i + 1.0, j + 0.5)).norm(),
                        (quad.vs<-1, +1>() - get_vertex(i + 0.0, j + 1.0)).norm(),
                        (quad.vs< 0, +1>() - get_vertex(i + 0.5, j + 1.0)).norm(),
                        (quad.vs<+1, +1>() - get_vertex(i + 1.0, j + 1.0)).norm(),

                });
    }

    // Линейный размер
    double H = std::max(hx, hy);

    cell_volume_error /= volume;
    cell_center_error /= H;

    face_area_error   /= H;
    face_center_error /= H;
    vertices_error    /= H;

    std::cout << "Errors\n";
    std::cout << "  Cell volume:  " << cell_volume_error << "\n";
    std::cout << "  Cell center:  " << cell_volume_error << "\n";
    std::cout << "  Face areas:   " << face_area_error << "\n";
    std::cout << "  Face centers: " << face_center_error << "\n";
    std::cout << "  Face normals: " << face_normal_error << "\n";
    std::cout << "  Vertices:     " << vertices_error << "\n";

    // Исправить сетку, если необходимо
    if (fix) {
        for (auto& cell: local_cells_) {
            index_t ic = cell.index();

            int i = std::round((cell.x() - 0.5 * hx - xmin) / hx);
            int j = std::round((cell.y() - 0.5 * hy - ymin) / hy);

            local_cells_.volume[ic] = volume;
            local_cells_.center[ic] = get_center(i, j);

            index_t iface = local_cells_.faces.offsets[ic];

            local_cells_.faces.area[iface + Side2D::L] = hy;
            local_cells_.faces.area[iface + Side2D::R] = hy;
            local_cells_.faces.area[iface + Side2D::B] = hx;
            local_cells_.faces.area[iface + Side2D::T] = hx;

            local_cells_.faces.normal[iface + Side2D::L] = -Vector3d::UnitX();
            local_cells_.faces.normal[iface + Side2D::R] =  Vector3d::UnitX();
            local_cells_.faces.normal[iface + Side2D::B] = -Vector3d::UnitY();
            local_cells_.faces.normal[iface + Side2D::T] =  Vector3d::UnitY();

            local_cells_.faces.center[iface + Side2D::L] = lface_center(i, j);
            local_cells_.faces.center[iface + Side2D::R] = rface_center(i, j);
            local_cells_.faces.center[iface + Side2D::B] = bface_center(i, j);
            local_cells_.faces.center[iface + Side2D::T] = tface_center(i, j);

            SqQuad& quad = local_cells_.verts.mapping<2>(ic);
            quad.vs<-1, -1>() = get_vertex(i + 0.0, j + 0.0);
            quad.vs< 0, -1>() = get_vertex(i + 0.5, j + 0.0);
            quad.vs<+1, -1>() = get_vertex(i + 1.0, j + 0.0);
            quad.vs<-1,  0>() = get_vertex(i + 0.0, j + 0.5);
            quad.vs< 0,  0>() = get_vertex(i + 0.5, j + 0.5);
            quad.vs<+1,  0>() = get_vertex(i + 1.0, j + 0.5);
            quad.vs<-1, +1>() = get_vertex(i + 0.0, j + 1.0);
            quad.vs< 0, +1>() = get_vertex(i + 0.5, j + 1.0);
            quad.vs<+1, +1>() = get_vertex(i + 1.0, j + 1.0);
        }
    }
}

std::string bytes(size_t n_bytes) {
    if (n_bytes < 1'000) {
        return std::format("{:6d} B ", n_bytes);
    }
    if (n_bytes < 1'000'000) {
        return std::format("{:6.1f} KB", 1.0e-3 * double(n_bytes));
    }
    if (n_bytes < 1'000'000'000) {
        return std::format("{:6.1f} MB", 1.0e-6 * double(n_bytes));
    }
    return std::format("{:6.1f} GB", 1.0e-9 * double(n_bytes));
}

void EuMesh::memory_usage() const {
    memory_t geom_size = local_cells_.memory_usage();
    memory_t face_size = local_cells_.faces.memory_usage();
    memory_t adj_size = local_cells_.faces.adjacent.memory_usage();
    memory_t vert_size = local_cells_.verts.memory_usage();
    memory_t node_size = local_nodes_.memory_usage();
    memory_t inc_size = local_nodes_.incident.memory_usage();

    memory_t cells_size;
    cells_size.needed = geom_size.needed + face_size.needed + adj_size.needed + vert_size.needed;
    cells_size.actual = geom_size.actual + face_size.actual + adj_size.actual + vert_size.actual;

    size_t geom_per_cell = cells_size.needed / local_cells_.n_cells();

    memory_t nodes_size;
    nodes_size.needed = node_size.needed + inc_size.needed;
    nodes_size.actual = node_size.actual + inc_size.actual;

    size_t geom_per_node = nodes_size.needed / local_nodes_.n_nodes();

    std::cout << "Local cells: " << local_cells_.n_cells() << "\n";
    std::cout << "  Cells:   " << bytes(cells_size.needed)  << " / " << bytes(cells_size.actual);
    std::cout << " (avg" << bytes(geom_per_cell) << " per cell)\n";
    std::cout << "    Geom:  " << bytes(geom_size.needed) << " / " << bytes(geom_size.actual) << "\n";
    std::cout << "    Adj:   " << bytes(adj_size.needed)  << " / " << bytes(adj_size.actual) << "\n";
    std::cout << "    Verts: " << bytes(vert_size.needed) << " / " << bytes(vert_size.actual) << "\n";
    std::cout << "    Faces: " << bytes(face_size.needed) << " / " << bytes(face_size.actual) << "\n";
    std::cout << "  Nodes:   " << bytes(nodes_size.needed)  << " / " << bytes(nodes_size.actual);
    std::cout << " (avg" << bytes(geom_per_node) << " per node)\n";
    std::cout << "    Geom:  " << bytes(node_size.needed) << " / " << bytes(node_size.actual) << "\n";
    std::cout << "    Inc:   " << bytes(inc_size.needed) << " / " << bytes(inc_size.actual) << "\n";

    memory_t data_size = local_cells_.data.memory_usage();
    size_t data_per_cell = data_size.needed / local_cells_.n_cells();
    std::cout << "  Data:    " << bytes(data_size.needed) << " / " << bytes(data_size.actual);
    std::cout << " (avg" << bytes(data_per_cell) << " per cell)\n";

    // TODO: Проверить размеры unique_nodes, tourists, migrants.
    // Сделать MPI-версию проверки памяти. Ну и подумать над оптимизацией размеров
}

int EuMesh::check_base() const {
    if (local_cells_.empty()) {
        if (mpi::single()) {
            std::cout << "\tEmpty storage\n";
            return -1;
        } else {
            return 0;
        }
    }

    auto dim = local_cells_.dim();

    if (dim != 2 && dim != 3) {
        std::cout << "\tDimension is not 2 or 3\n";
        return -1;
    }

    int res = 0;
    for (index_t ic = 0; ic < local_cells_.n_cells(); ++ic) {
        if (local_cells_.index[ic] < 0 || local_cells_.index[ic] != ic) {
            std::cout << "\tWrong cell index\n";
            return -1;
        }

        if (local_cells_.rank[ic] < 0 || local_cells_.rank[ic] != mpi::rank()) {
            std::cout << "\tWrong cell rank\n";
            return -1;
        }

        // Проверим число вершин
        int n_nodes = local_cells_.verts.count(ic);
        int n_max_nodes = local_cells_.verts.max_count(ic);
        if (local_cells_.adaptive()) {
            if ((dim == 2 && n_nodes == n_max_nodes && n_max_nodes != 9) ||
                (dim == 3 && n_nodes == n_max_nodes && n_max_nodes != 27)) {
                std::cout << "\tCell has wrong node count " << n_nodes << " " << n_max_nodes << "\n";
                local_cells_.print_info(ic);
                return -1;
            }
        }
        else {
            if (n_nodes != n_max_nodes) {
                std::cout << "\tCell has strange number of nodes (" << n_nodes << ")\n";
                local_cells_.print_info(ic);
                return -1;
            }
            if (n_nodes < dim + 1) {
                std::cout << "\tCell has too little nodes (" << n_nodes << ")\n";
                local_cells_.print_info(ic);
                return -1;
            }
        }

        // Проверим число граней
        int n_faces = local_cells_.face_count(ic);
        int n_max_faces = local_cells_.faces.max_count(ic);
        if (local_cells_.adaptive()) {
            for (int i = 0; i < FpC(dim); ++i) {
                if (local_cells_.faces.is_undefined(local_cells_.faces.offsets[ic] + i)) {
                    std::cout << "\tCell has no one of main faces\n";
                    local_cells_.print_info(ic);
                    return -1;
                }
            }
            if (n_faces > FpC(dim)) {
                std::cout << "\tCell has too much faces (" << local_cells_.face_count(ic) << ")\n";
                local_cells_.print_info(ic);
                return -1;
            }
            if ((dim == 2 && (n_faces > n_max_faces || n_max_faces != 8)) ||
                (dim == 3 && (n_faces > n_max_faces || n_max_faces != 24))) {
                std::cout << "\tCell has wrong face count " << n_faces << " " << n_max_faces << "\n";
                local_cells_.print_info(ic);
                return -1;
            }
        }
        else {
            if (n_faces != n_max_faces) {
                std::cout << "\tCell has strange number of faces (" << n_faces << ")\n";
                local_cells_.print_info(ic);
                return -1;
            }
            if (n_faces < dim + 1) {
                std::cout << "\tCell has too little faces (" << n_max_faces << ")\n";
                local_cells_.print_info(ic);
                return -1;
            }
            for (int i = 0; i < n_faces; ++i) {
                if (local_cells_.faces.is_undefined(local_cells_.faces.offsets[ic] + i)) {
                    std::cout << "\tCell has undefined face\n";
                    local_cells_.print_info(ic);
                    return -1;
                }
            }
        }

        // Правильное задание геометрии
        res = local_cells_.check_geometry(ic);
        if (res < 0) return res;

        // Грани правильно ориентированы
        res = local_cells_.check_base_face_orientation(ic);
        if (res < 0) return res;

        // Порядок основных вершин
        res = local_cells_.check_base_vertices_order(ic);
        if (res < 0) return res;

        // Проверка смежности
#ifndef ZEPHYR_MPI
        res = local_cells_.check_connectivity(ic);
#else
        res = local_cells_.check_connectivity(ic, tourists_.ghost_cells());
#endif
        if (res < 0) return res;
    }

    // Если уникальные узлы не построены, то завершаем
    if (!local_cells_.verts.has_nodes()) {
        if (!local_cells_.verts.index.empty() || !local_cells_.verts.ghost.empty()) {
            std::cout << "\tUnique nodes: not empty verts arrays\n";
            return -1;
        }
        return 0;
    }

#ifndef ZEPHYR_MPI
    res = local_nodes_.check_nodes(local_cells_);
#else
    res = local_nodes_.check_nodes(local_cells_,
        tourists_.ghost_cells(), tourists_.ghost_nodes());
#endif
    if (res < 0) return res;

    return 0;
}

int EuMesh::check_refined() const {
    if (local_cells_.empty()) {
        if (mpi::single()) {
            std::cout << "\tEmpty storage\n";
            return -1;
        } else {
            return 0;
        }
    }

    auto dim = local_cells_.dim();

    if (dim != 2 && dim != 3) {
        std::cout << "\tDimension is not 2 or 3\n";
        return -1;
    }

    int res = 0;
    for (index_t ic = 0; ic < local_cells_.n_cells(); ++ic) {
        if (local_cells_.is_undefined(ic)) {
            std::cout << "\tUndefined cell\n";
            return -1;
        }

        if (local_cells_.index[ic] < 0 || local_cells_.index[ic] != ic) {
            std::cout << "\tWrong cell index\n";
            return -1;
        }

        if (local_cells_.rank[ic] < 0 || local_cells_.rank[ic] != mpi::rank()) {
            std::cout << "\tWrong cell rank\n";
            return -1;
        }

        // Число граней
        for (int i = 0; i < FpC(dim); ++i) {
            if (local_cells_.faces.is_undefined(local_cells_.faces.offsets[ic] + i)) {
                std::cout << "\tCell has no one of main faces\n";
                local_cells_.print_info(ic);
                return -1;
            }
        }

        // Вершины дублируются
        for (int i = local_cells_.verts.offsets[ic]; i < local_cells_.verts.offsets[ic + 1]; ++i) {
            for (int j = i + 1; j < local_cells_.verts.offsets[ic + 1]; ++j) {
                double dist = (local_cells_.verts[i] - local_cells_.verts[j]).norm();
                if (dist < 1.0e-5 * local_cells_.linear_size(ic)) {
                    std::cout << "\tIdentical vertices\n";
                    local_cells_.print_info(ic);
                    return -1;
                }
            }
        }

        // Правильное задание геометрии
        res = local_cells_.check_geometry(ic);
        if (res < 0) return res;

        // Грани правильно ориентированы
        res = local_cells_.check_base_face_orientation(ic);
        if (res < 0) return res;

        // Порядок основных вершин
        res = local_cells_.check_base_vertices_order(ic);
        if (res < 0) return res;

        // Проверка сложных граней
        res = local_cells_.check_complex_faces(ic);
        if (res < 0) return res;

        // Проверка смежности
#ifndef ZEPHYR_MPI
        res = local_cells_.check_connectivity(ic);
#else
        res = local_cells_.check_connectivity(ic, tourists_.ghost_cells());
#endif
        if (res < 0) return res;
    }

    // Если уникальные узлы не построены, то завершаем
    if (!local_cells_.verts.has_nodes()) {
        if (!local_cells_.verts.index.empty() || !local_cells_.verts.ghost.empty()) {
            std::cout << "\tUnique nodes: not empty verts arrays\n";
            return -1;
        }
        return 0;
    }

#ifndef ZEPHYR_MPI
    res = local_nodes_.check_nodes(local_cells_);
#else
    res = local_nodes_.check_nodes(local_cells_,
        tourists_.ghost_cells(), tourists_.ghost_nodes());
#endif
    if (res < 0) return res;

    return 0;
}

Box EuMesh::bbox() const {
    Box box1 = Box::Empty(3);
    for (auto& v: local_cells_.verts.coord) {
        box1.capture(v);
    }

    Box box2(box1);

#ifdef ZEPHYR_MPI
    if (!mpi::single()) {
        // Покомпонентный минимум/максимум
        MPI_Allreduce(box1.vmin.data(), box2.vmin.data(), 3, MPI_DOUBLE, MPI_MIN, mpi::comm());
        MPI_Allreduce(box1.vmax.data(), box2.vmax.data(), 3, MPI_DOUBLE, MPI_MAX, mpi::comm());
    }
#endif

    return box2;
}

void EuMesh::push_back(const geom::Line& line) {
    geom::Polygon poly = {line[0], line[1], line[1], line[0]};
    local_cells_.push_back(poly);
}

void EuMesh::push_back(const geom::Polygon& poly) {
    local_cells_.push_back(poly);
}

void EuMesh::push_back(const geom::Polyhedron& poly) {
    local_cells_.push_back(poly);
}

void EuMesh::add_marker(const geom::Vector3d& pos, double size) {
    if (dim() < 3) {
        double c1 = 0.5 * size;
        double c2 = 0.5 * size * std::sqrt(3.0);
        Vector3d a = {0.0, -size, 1.0};
        Vector3d b = {+c2, c1, 1.0};
        Vector3d c = {-c2, c1, 1.0};
        Polygon poly = {pos + a, pos + b, pos + c};
        push_back(poly);
    }
    else {
        throw std::runtime_error("3D markers not supported");
    }
}

EuCell_Iter EuMesh::begin() {
    return {&local_cells_, 0,
        mpi_cond(&tourists_.ghost_cells(), nullptr) };
}

EuCell_Iter EuMesh::end() {
    return {&local_cells_, local_cells_.n_cells(),
        mpi_cond(&tourists_.ghost_cells(), nullptr) };
}

EuCell EuMesh::operator[](index_t idx) {
    return {&local_cells_, idx,
        mpi_cond(&tourists_.ghost_cells(), nullptr) };
}

EuCell EuMesh::operator()(int i, int j) {
    i = (i + nx_) % nx_;
    j = (j + ny_) % ny_;
    return operator[](ny_ * i + j);
}

EuCell EuMesh::operator()(int i, int j, int k) {
    i = (i + nx_) % nx_;
    j = (j + ny_) % ny_;
    k = (k + nz_) % nz_;
    return operator[](nz_ * (ny_ * i + j) + k);
}

EuNodeRange EuMesh::nodes() {
    if (!has_nodes()) {
        throw std::runtime_error("EuMesh::nodes: has no unique nodes");
    }
#ifndef ZEPHYR_MPI
    return EuNodeRange(&local_nodes_, &local_cells_, nullptr);
#else
    return EuNodeRange(&local_nodes_, &local_cells_, &tourists_.ghost_cells());
#endif
}

void EuMesh::backup(const std::string& sroot, const std::vector<std::string>& variables) const {
    namespace fs = std::filesystem;

    const fs::path root = sroot;

    if (mpi::master()) {
        if (fs::exists(root)) {
            fs::remove_all(root);
        }
        fs::create_directories(root);
    }

    // Основной json
    std::ofstream file;
    if (mpi::master()) {
        file.open(root / "backup.json", std::ios::out | std::ios::trunc);
        file << "{\n";
        file << "  \"mesh\": {\n";
    }

    if (mpi::master()) {
        file << std::boolalpha;
        file << "    \"max_level\":  " << max_level_ << ",\n";
        file << "    \"structured\": " << structured_ << ",\n";
        if (structured_) {
            file << "    \"nx\": " << nx_ << ",\n";
            file << "    \"ny\": " << ny_ << ",\n";
            file << "    \"nz\": " << nz_ << ",\n";
        }
    }

    if (mpi::master()) {
        fs::create_directory(root / "cells");
        file << "    \"cells\": {\n";
    }

    local_cells_.backup(root, file, "      ", variables);

    if (mpi::master()) {
        file << "    }\n"; // mesh.cells
        file << "  }\n"; // mesh
        file << "}";
        file.close();
    }
}

} // namespace zephyr::mesh