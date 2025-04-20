#include <fstream>

#include <zephyr/mesh/euler/soa_mesh.h>

#include <zephyr/mesh/amr2/apply.h>
#include <zephyr/mesh/amr2/balancing.h>

namespace zephyr::mesh {

void SoaCell::resize(index_t n_cells,
    index_t faces_per_cell, index_t nodes_per_cell) {
    // Поля ячеек по числу ячеек, логично
    rank.resize(n_cells, -1);
    index.resize(n_cells, -1);
    next.resize(n_cells, -1);

    center.resize(n_cells);
    volume.resize(n_cells);
    volume_alt.resize(n_cells);

    // +1 для заключительной
    face_begin.resize(n_cells + 1);
    node_begin.resize(n_cells + 1);

    flag.resize(n_cells);
    b_idx.resize(n_cells);
    z_idx.resize(n_cells);
    level.resize(n_cells);

    // Поля граней и вершин пока так
    faces.resize(n_cells * faces_per_cell);
    verts.resize(n_cells * nodes_per_cell);

    // Поля данных только для ячеек
    data.resize(n_cells);
}

inline double trimin(double x, double y, double z) {
    return std::min(x, std::min(y, z));
}

double QCell::diameter() const {
    if (m_cells->adaptive) {
        if (m_cells->dim == 2) {
            const SqQuad& vertices = m_cells->get_vertices<2>(cell_idx);
            return std::sqrt(std::min(
                    (vertices.vs<+1, 0>() - vertices.vs<-1, 0>()).squaredNorm(),
                    (vertices.vs<0, +1>() - vertices.vs<0, -1>()).squaredNorm()));
        } else {
            const SqCube& vertices = m_cells->get_vertices<3>(cell_idx);
            return std::sqrt(trimin(
                    (vertices.vs<+1, 0, 0>() - vertices.vs<-1, 0, 0>()).squaredNorm(),
                    (vertices.vs<0, +1, 0>() - vertices.vs<0, -1, 0>()).squaredNorm(),
                    (vertices.vs<0, 0, +1>() - vertices.vs<0, 0, -1>()).squaredNorm()));
        }
    }
    else {
        throw std::runtime_error("Not implemented");
    }
}

bool FaceIt::to_skip(Direction dir) const {
    if (m_cells->faces.boundary[iface] == Boundary::UNDEFINED) {
        return true;
    }
    switch (dir) {
        case Direction::ANY:
            return false;
        case Direction::X:
            return std::abs(m_cells->faces.normal[iface].x()) < 0.7;
        case Direction::Y:
            return std::abs(m_cells->faces.normal[iface].y()) < 0.7;
        case Direction::Z:
            return std::abs(m_cells->faces.normal[iface].z()) < 0.7;
        default:
            return false;
    }
}

FacesIts QCell::faces(Direction dir) {
    return FacesIts(m_cells,
                    m_cells->face_begin[cell_idx],
                    m_cells->face_begin[cell_idx + 1],
                    dir);
}

bool QFace::is_boundary() const {
    return m_cells->faces.boundary[iface] != Boundary::ORDINARY &&
           m_cells->faces.boundary[iface] != Boundary::PERIODIC &&
           m_cells->faces.boundary[iface] != Boundary::UNDEFINED;
}

bool QFace::is_actual() const {
    return m_cells->faces.boundary[iface] != Boundary::UNDEFINED;
}

bool QFace::is_undefined() const {
    return  m_cells->faces.boundary[iface] == Boundary::UNDEFINED;
}

void QFace::set_undefined() {
    m_cells->faces.boundary[iface] = Boundary::UNDEFINED;
    m_cells->faces.adjacent[iface].rank  = -1;
    m_cells->faces.adjacent[iface].index = -1;
    m_cells->faces.adjacent[iface].alien = -1;
}

Vector3d QFace::symm_point(const Vector3d& p) const {
    return p + 2.0 * (m_cells->faces.center[iface] - p).dot(m_cells->faces.normal[iface]) * m_cells->faces.normal[iface];
}

double QFace::area() const { return m_cells->faces.area[iface]; }

double QFace::area(bool axial) const {
    return axial ? m_cells->faces.area_alt[iface] : m_cells->faces.area[iface];
}

QCell QFace::neib() const {
    return {m_cells, neib_index()};
}

index_t QFace::neib_index() const { return m_cells->faces.adjacent[iface].index; }

Vector3d QFace::neib_center() const { return m_cells->center[neib_index()]; }

Boundary QFace::flag() const { return m_cells->faces.boundary[iface]; }

const Vector3d& QFace::normal() const { return m_cells->faces.normal[iface]; }

const Vector3d& QFace::center() const { return m_cells->faces.center[iface]; }

double QFace::x() const { return center().x(); }

inline double QFace::y() const { return center().y(); }

inline double QFace::z() const { return center().z(); }

inline const Adjacent &QFace::adjacent() const {
    return m_cells->faces.adjacent[iface];
}


SoaMesh::SoaMesh(EuMesh &mesh)
    : cells(mesh.locals()) {

}

SoaMesh::SoaMesh(Generator &gen) {
    EuMesh mesh(gen, int{});
    cells.initialize(mesh.locals());
}

SoaCell::SoaCell(AmrStorage &locals) {
    initialize(locals);
}

void SoaCell::initialize(AmrStorage &locals) {
    n_cells = locals.size();

    if (locals.empty()) {
        throw std::runtime_error("SoaMesh: locals is empty");
    }

    dim = locals[0].dim;
    adaptive = locals[0].adaptive;
    linear = locals[0].linear;
    axial = locals[0].axial;

    index_t faces_per_cell = dim < 3 ? 8 : 24;
    index_t nodes_per_cell = dim < 3 ? 9 : 27;

    resize(n_cells, faces_per_cell, nodes_per_cell);


    for (index_t ic = 0; ic < n_cells; ++ic) {
        rank[ic] = locals[ic].rank;
        index[ic] = locals[ic].index;
        next[ic] = locals[ic].next;

        center[ic] = locals[ic].center;
        volume[ic] = locals[ic].volume;
        volume_alt[ic] = locals[ic].volume_alt;
        face_begin[ic] = ic * faces_per_cell;
        node_begin[ic] = ic * nodes_per_cell;

        for (index_t jf = 0; jf < faces_per_cell; ++jf) {
            index_t iface = ic * faces_per_cell + jf;

            faces.boundary[iface] = locals[ic].faces[jf].boundary;
            faces.adjacent[iface] = locals[ic].faces[jf].adjacent;
            faces.normal[iface] = locals[ic].faces[jf].normal;
            faces.center[iface] = locals[ic].faces[jf].center;
            faces.area[iface] = locals[ic].faces[jf].area;
            faces.area_alt[iface] = locals[ic].faces[jf].area_alt;
            faces.vertices[iface] = locals[ic].faces[jf].vertices;
        }

        for (index_t jn = 0; jn < nodes_per_cell; ++jn) {
            index_t v_idx = ic * nodes_per_cell + jn;
            verts[v_idx] = locals[ic].vertices[jn];
        }
    }
    face_begin[n_cells] = n_cells * faces_per_cell;
    node_begin[n_cells] = n_cells * nodes_per_cell;

    data.resize(n_cells);
}

void SoaMesh::init_amr() {
    // Эта функция нифига не используется же?
    m_max_level = 0;

    if (cells.empty()) {
        return;
    }

    int rank = mpi::rank();
    int size = mpi::size();

    auto cells_nums = mpi::all_gather(cells.size());
    std::vector<int> offset(size, 0);
    for (int r = 0; r < size - 1; ++r) {
        offset[r + 1] = offset[r] + cells_nums[r];
    }

    for (int ic = 0; ic < cells.size(); ++ic) {
        cells.b_idx[ic] = offset[rank] + ic;
        cells.z_idx[ic] = 0;
        cells.level[ic] = 0;
        cells.flag[ic] = 0;
    }
}

bool SoaMesh::is_adaptive() const {
    return m_max_level > 0;
}

int SoaMesh::max_level() const {
    return m_max_level;
}

void SoaMesh::set_max_level(int max_level) {
    m_max_level = std::max(0, std::min(max_level, 15));
}

void SoaMesh::set_distributor(Distributor distr) {
    distributor = std::move(distr);
}

void SoaMesh::balance_flags() {
#if SCRUTINY
    static bool first_time = true;
    if (first_time) {
        int res = check_base();
        if (res < 0) {
            throw std::runtime_error("Check base failed");
        }
        first_time = false;
    }
#endif
    if (mpi::single()) {
        std::cerr << "No balancing sorry\n";
        //amr::balance_flags(cells, m_max_level);
    }
#ifdef ZEPHYR_MPI
    else {
        amr::balance_flags(m_locals, m_aliens, m_max_level, *this);
    }
#endif
}

void SoaMesh::apply_flags() {
    if (mpi::single()) {
        std::cerr << "No apply sorry\n";
        //amr::apply(cells, distributor);
    }
#ifdef ZEPHYR_MPI
    else {
        amr::apply(m_locals, m_aliens, distributor, *this);
    }
#endif

#if SCRUTINY
    int res = check_refined();
    if (res < 0) {
        throw std::runtime_error("Check refined failed");
    }
#endif
}

void SoaMesh::refine() {
    if (!is_adaptive()) { return; }

    static Stopwatch balance;
    static Stopwatch apply;
    static Stopwatch full;

    // Для однопроцессорной версии при пустой сетке сразу выход
    if (mpi::single() && cells.empty()) {
        throw std::runtime_error("EuMesh::refine() error: Empty mesh");
    }

    full.resume();

    balance.resume();
    balance_flags();
    balance.stop();

    apply.resume();
    apply_flags();
    apply.stop();

    full.stop();

#if CHECK_PERFORMANCE
    static size_t counter = 0;
    if (counter % amr::check_frequency == 0) {
        std::cout << "  Balance steps elapsed: " << std::setw(10) << balance.milliseconds() << " ms\n";
        std::cout << "  Apply steps elapsed:   " << std::setw(10) << apply.milliseconds() << " ms\n";
        std::cout << "Refine steps elapsed:    " << std::setw(10) << full.milliseconds() << " ms\n";
    }
    ++counter;
#endif
}


void SoaCell::print_info(index_t ic) const {
    std::cout << "\t\tcenter: " << center[ic].transpose() << "\n";
    std::cout << "\t\trank:   " << rank[ic] << "\n";
    std::cout << "\t\tindex:  " << index[ic] << "\n";
    std::cout << "\t\tsize:   " << linear_size(ic) << "\n";
    std::cout << "\t\tvolume: " << volume[ic] << "\n";
    std::cout << "\t\tflag:   " << flag[ic] << "\n";
    std::cout << "\t\tnext:   " << next[ic] << "\n";
    std::cout << "\t\tb_idx:  " << b_idx[ic] << "\n";
    std::cout << "\t\tlevel:  " << level[ic] << "\n";
    std::cout << "\t\tz_idx:  " << z_idx[ic] << "\n";

    std::cout << "\t\tcell.vertices:\n";
    for (int i = node_begin[ic]; i < node_begin[ic + 1]; ++i) {
        if (!verts[i].hasNaN()) {
            std::cout << "\t\t\t" << i << ": " << verts[i].transpose() << "\n";
        }
    }

    std::cout << "\t\tcell.faces:\n";
    for (int iface = face_begin[ic]; iface < face_begin[ic + 1]; ++iface) {
        if (faces.is_undefined(iface)) continue;

        std::cout << "\t\t\t" << side_to_string(iface) << ":\n";
        std::cout << "\t\t\t\tvertices:";
        for (int j = 0; j < (dim < 3 ? 2 : 4); ++j) {
            std::cout << " " << faces.vertices[iface][j];
        }
        std::cout << "\n";
        std::cout << "\t\t\t\tflag:       " << faces.boundary[iface] << "\n";
        std::cout << "\t\t\t\tarea:       " << faces.area[iface] << "\n";
        std::cout << "\t\t\t\tnormal:     " << faces.normal[iface].transpose() << "\n";
        std::cout << "\t\t\t\tadj.rank:   " << faces.adjacent[iface].rank << "\n";
        std::cout << "\t\t\t\tadj.index:  " << faces.adjacent[iface].index << "\n";
        std::cout << "\t\t\t\tadj.alien:  " << faces.adjacent[iface].alien << "\n";
    }
}

void SoaCell::visualize(index_t ic, std::string filename) const {
    if (filename.find(".py") == std::string::npos) {
        filename += ".py";
    }

    std::ofstream file(filename);

    file << std::scientific << std::setprecision(6);

    file << "#!/bin/python3\n";
    file << "import numpy as np\n";
    file << "import matplotlib.pyplot as plt\n\n\n";

    file << "def spline(v1, vc, v2, t):\n";
    file << "    return vc + 0.5*(v2 - v1)*t + 0.5*(v2 - 2*vc + v1)*t**2\n\n\n";

    file << "def plot_arrow(ax, x, y, vx, vy, color, width=0.5):\n";
    file << "    phi = 0.1*np.pi\n";
    file << "    ax.plot([x - (+vx*np.cos(phi) + vy*np.sin(phi)), x, x - (vx*np.cos(phi) - vy*np.sin(phi))],\n"
         << "            [y - (-vx*np.sin(phi) + vy*np.cos(phi)), y, y - (vx*np.sin(phi) + vy*np.cos(phi))],\n"
         << "            color=color)\n\n";

    // Основные точки
    SqQuad map;
    for (int i = node_begin[ic]; i < node_begin[ic + 1]; ++i) {
        map[i] = verts[i];
    }

    file << "fig = plt.figure(dpi=150, figsize=(8, 8))\n";
    file << "ax = fig.add_subplot()\n\n";
    file << "ax.set_aspect('equal')\n";

    file << "ax.plot([";
    for (int i = 0; i < 8; ++i) {
        file << map[i].x() << ", ";
    }
    file << map[8].x() << "],\n";
    file << "        [";
    for (int i = 0; i < 8; ++i) {
        file << map[i].y() << ", ";
    }
    file << map[8].y() << "],\n";
    file << "        linestyle='none', color='orange', marker='o')\n\n";

    file << "ax.plot([" << center[ic].x() << "], [" << center[ic].y() << "], color='black', marker='x')\n\n";

    if (dim != 2) {
        throw std::runtime_error("Can't visualize 3D cell, sorry");
    }

    SqQuad vertices = get_vertices<2>(ic);

    for (int i = 0; i < 9; ++i) {
        file << "ax.text(" << vertices[i].x() << ", " << vertices[i].y() << ", " << i << ")\n";
    }

    file << "unit = np.linspace(-1.0, 1.0, 101)\n\n";

    SqLine L(SqLine({vertices[0], vertices[3], vertices[6]}));
    SqLine R(SqLine({vertices[2], vertices[5], vertices[8]}));
    SqLine B(SqLine({vertices[0], vertices[1], vertices[2]}));
    SqLine T(SqLine({vertices[6], vertices[7], vertices[8]}));

    for (auto& sf: {L, R, B, T}) {
        file << "curve_Lx = spline(" << sf(-1).x() << ", " << sf(0).x() << ", " << sf(+1).x() << ", unit)\n";
        file << "curve_Ly = spline(" << sf(-1).y() << ", " << sf(0).y() << ", " << sf(+1).y() << ", unit)\n";
        file << "ax.plot(curve_Lx, curve_Ly, linestyle='dotted', color='green', linewidth=0.5)\n\n";
    }

    for (int iface = face_begin[ic]; iface < face_begin[ic + 1]; ++iface) {
        if (faces.is_undefined(iface)) {
            continue;
        }

        double area = faces.area[iface];
        Vector3d normal = faces.normal[iface];

        Vector3d v1 = vertices[faces.vertices[iface][0]];
        Vector3d v2 = vertices[faces.vertices[iface][1]];
        Vector3d vc = faces.center[iface];

        file << "# face " << iface << "\n";

        file << "ax.plot([" << v1.x() << ", " << v2.x() << "], ["
             << v1.y() << ", " << v2.y() << "], color='green', marker='.')\n";

        double a = 0.03;
        file << "plot_arrow(ax, " << vc.x() << ", " << vc.y() << ", "
             << a * (v2.x() - v1.x()) << ", " << a * (v2.y() - v1.y()) << ", color='green')\n";

        double b = 0.25 * area;
        file << "ax.plot([" << vc.x() << ", " << vc.x() + b * normal.x() << "], ["
             << vc.y() << ", " << vc.y() + b * normal.y() << "], color='red')\n";

        double c = 0.2 * b;
        file << "plot_arrow(ax, " << vc.x() + b * normal.x() << ", " << vc.y() + b * normal.y() << ", "
             << c * normal.x() << ", " << c * normal.y() << ", color='red')\n\n";
    }

    file << "\n";
    file << "fig.tight_layout()\n";
    file << "plt.show()\n";
}

int SoaCell::check_geometry(index_t ic) const {
    for (index_t iface = face_begin[ic]; iface < face_begin[ic + 1]; ++iface) {
        if (faces.is_undefined(iface)) continue;

        Vector3d fc(0.0, 0.0, 0.0);
        for (int iv = 0; iv < VpF(dim); ++iv) {
            fc += verts[node_begin[ic] + faces.vertices[iface][iv]];
        }
        fc /= VpF(dim);

        // Нормаль внешняя
        if (faces.normal[iface].dot(fc - center[ic]) < 0.0) {
            std::cout << "\tWrong normal direction (inside cell)\n";
            print_info(ic);
            return -1;
        }

        // Вершины грани перечислены в правильном порядке
        if (dim > 2) {
            Vector3d v0 = verts[node_begin[ic] + faces.vertices[iface][0]];
            Vector3d v1 = verts[node_begin[ic] + faces.vertices[iface][1]];
            Vector3d v2 = verts[node_begin[ic] + faces.vertices[iface][2]];
            Vector3d v3 = verts[node_begin[ic] + faces.vertices[iface][3]];

            Vector3d n1 = (v2 - v1).cross(v0 - v1);
            Vector3d n2 = (v1 - v2).cross(v3 - v2);
            if (n1.dot(n2) < 0.0) {
                std::cout << "Wrong order of vertices on face\n";
                print_info(ic);
                return -1;
            }
        }
    }
    return 0;
}

int SoaCell::check_base_face_orientation(index_t ic) const {
    if (dim == 2) {
        Vector3d nx1 = faces.normal[face_begin[ic] + Side::L];
        Vector3d nx2 = faces.normal[face_begin[ic] + Side::R];
        Vector3d ny1 = faces.normal[face_begin[ic] + Side::B];
        Vector3d ny2 = faces.normal[face_begin[ic] + Side::T];

        if (nx1.dot(nx2) > -0.8) {
            std::cout << "\tOpposite outward normals (left-right) are co-directed\n";
            print_info(ic);
            return -1;
        }
        if (ny1.dot(ny2) > -0.8) {
            std::cout << "\tOpposite outward normals (bottom-top) are co-directed\n";
            print_info(ic);
            return -1;
        }
        if (nx1.cross(ny1).z() < 0.0) {
            std::cout << "\tWrong face orientation (left-bottom)\n";
            print_info(ic);
            return -1;
        }
        if (nx2.cross(ny2).z() < 0.0) {
            std::cout << "\tWrong face orientation (right-top)\n";
            print_info(ic);
            return -1;
        }
    } else {
        Vector3d nx1 = faces.normal[face_begin[ic] + Side::L];
        Vector3d nx2 = faces.normal[face_begin[ic] + Side::R];
        Vector3d ny1 = faces.normal[face_begin[ic] + Side::B];
        Vector3d ny2 = faces.normal[face_begin[ic] + Side::T];
        Vector3d nz1 = faces.normal[face_begin[ic] + Side::X];
        Vector3d nz2 = faces.normal[face_begin[ic] + Side::F];

        if (nx1.dot(nx2) > -0.8) {
            std::cout << "\tOpposite outward normals (left-right) are co-directed\n";
            print_info(ic);
            return -1;
        }
        if (ny1.dot(ny2) > -0.8) {
            std::cout << "\tOpposite outward normals (bottom-top) are co-directed\n";
            print_info(ic);
            return -1;
        }
        if (nz1.dot(nz2) > -0.8) {
            std::cout << "\tOpposite outward normals (back-front) are co-directed\n";
            print_info(ic);
            return -1;
        }
        if (nx1.dot(ny1.cross(nz1)) > 0.0) {
            std::cout << "\tWrong face orientation (left-bottom-back)\n";
            print_info(ic);
            return -1;
        }
        if (nx2.dot(ny2.cross(nz2)) < 0.0) {
            std::cout << "\tWrong face orientation (right-top-front)\n";
            print_info(ic);
            return -1;
        }
    }
    return 0;
}

int SoaCell::check_base_vertices_order(index_t ic) const {
    if (dim == 2) {
        SqQuad verts = get_vertices<2>(ic);

        bool bad = false;

        // Индекс пересечения двух граней
        auto cross_face = [this](index_t iface_base, Side side1, Side side2) -> int {
            auto iface = iface_base + side1;
            auto jface = iface_base + side2;
            if (faces.is_undefined(iface)) {
                return 100;
            }
            if (faces.is_undefined(jface)) {
                return 100;
            }
            for (int i: {0, 1}) {
                for (int j: {0, 1}) {
                    if (faces.vertices[iface][i] == faces.vertices[jface][j]) {
                        return faces.vertices[iface][i];
                    }
                }
            }
            return 100;
        };

        for (int i: {-1, 0}) {
            for (int j: {-1, 0}) {
                Vector3d a = verts(i + 1, j) - verts(i, j);
                Vector3d b = verts(i, j + 1) - verts(i, j);

                Vector3d c = verts(i, j + 1) - verts(i + 1, j + 1);
                Vector3d d = verts(i + 1, j) - verts(i + 1, j + 1);

                if (a.cross(b).z() < 0.0 ||
                    c.cross(d).z() < 0.0) {
                    bad = true;
                    break;
                }
            }
            if (bad) {
                break;
            }
        }

        // Пересечения граней по нужным вершинам
        if (cross_face(face_begin[ic], Side::LEFT0, Side::BOTTOM0) != SqQuad::iss<-1, -1>()) {
            bad = true;
        }
        if (cross_face(face_begin[ic], Side::LEFT0, Side::TOP) != SqQuad::iss<-1, +1>() &&
            cross_face(face_begin[ic], Side::LEFT1, Side::TOP) != SqQuad::iss<-1, +1>()) {
            bad = true;
        }
        if (cross_face(face_begin[ic], Side::RIGHT, Side::BOTTOM0) != SqQuad::iss<+1, -1>() &&
            cross_face(face_begin[ic], Side::RIGHT, Side::BOTTOM1) != SqQuad::iss<+1, -1>()) {
            bad = true;
        }
        if (cross_face(face_begin[ic], Side::RIGHT0, Side::TOP0) != SqQuad::iss<+1, +1>() &&
            cross_face(face_begin[ic], Side::RIGHT0, Side::TOP1) != SqQuad::iss<+1, +1>() &&
            cross_face(face_begin[ic], Side::RIGHT1, Side::TOP0) != SqQuad::iss<+1, +1>() &&
            cross_face(face_begin[ic], Side::RIGHT1, Side::TOP1) != SqQuad::iss<+1, +1>()) {
            bad = true;
        }

        if (bad) {
            std::cout << "\tBad arrangement of vertices in cell\n";
            print_info(ic);
            return -1;
        }
    } else {
        SqCube verts = get_vertices<3>(ic);

        bool bad = false;

        // Индекс пересечения трех граней
        auto cross_face = [](const BFaces& faces, Side side1, Side side2, Side side3) -> int {
            auto& face1 = faces[side1];
            auto& face2 = faces[side2];
            auto& face3 = faces[side3];

            if (face1.is_undefined()) {
                return 100;
            }
            if (face2.is_undefined()) {
                return 100;
            }
            if (face3.is_undefined()) {
                return 100;
            }
            for (int i: {0, 1, 2, 3}) {
                for (int j: {0, 1, 2, 3}) {
                    for (int k: {0, 1, 2, 3}) {
                        if (face1.vertices[i] == face2.vertices[j] &&
                            face2.vertices[j] == face3.vertices[k]) {
                            return face1.vertices[i];
                        }
                    }
                }
            }
            return 100;
        };

        for (int i: {-1, 0}) {
            for (int j: {-1, 0}) {
                for (int k: {-1, 0}) {
                    Vector3d a = verts(i + 1, j, k) - verts(i, j, k);
                    Vector3d b = verts(i, j + 1, k) - verts(i, j, k);
                    Vector3d c = verts(i, j, k + 1) - verts(i, j, k);

                    Vector3d A = verts(i + 1, j + 1, k + 1) - verts(i, j + 1, k + 1);
                    Vector3d B = verts(i + 1, j + 1, k + 1) - verts(i + 1, j, k + 1);
                    Vector3d C = verts(i + 1, j + 1, k + 1) - verts(i + 1, j + 1, k);

                    if (a.dot(b.cross(c)) < 0.0 || A.dot(B.cross(C)) < 0.0) {
                        bad = true;
                        break;
                    }
                }
                if (bad) {
                    break;
                }
            }
            if (bad) {
                break;
            }
        }

        // Пересечения граней по нужным вершинам ???

        if (bad) {
            std::cout << "\tBad arrangement of vertices in cell\n";
            print_info(ic);
            return -1;
        }
    }

    return 0;
}

int SoaCell::check_complex_faces(index_t ic) const {
    if (dim == 2) {
        for (int s = 0; s < 4; ++s) {
            auto iface1 = face_begin[ic] + s;
            auto iface2 = face_begin[ic] + s + 6;
            if (faces.is_undefined(iface2)) {
                continue;
            }

            if (std::abs(faces.normal[iface1].dot(faces.normal[iface2]) - 1.0) > 1.0e-2) {
                std::cout << "\tSubfaces are not co-directed (" + side_to_string(s) + " side)\n";
                print_info(ic);
                return -1;
            }
            int count = 0;
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    if (faces.vertices[iface1][i] == faces.vertices[iface2][j]) {
                        ++count;
                    }
                }
            }
            if (count != 1) {
                std::cout << "\tStrange complex face\n";
                print_info(ic);
                return -1;
            }
        }
    } else {
        for (int s = 0; s < 6; ++s) {
            auto iface1 = face_begin[ic] + s;
            auto iface2 = face_begin[ic] + s + 6;
            if (faces.is_undefined(iface2)) {
                continue;
            }

            auto iface3 = face_begin[ic] + s + 12;
            auto iface4 = face_begin[ic] + s + 18;
            if (faces.is_undefined(iface3) || faces.is_undefined(iface4)) {
                std::cout << "\tComplex 3D face (" + side_to_string(s) + " side) has less than 4 subfaces\n";
                print_info(ic);
                return -1;
            }

            double d1 = std::abs(faces.normal[iface1].dot(faces.normal[iface2]) - 1.0);
            double d2 = std::abs(faces.normal[iface1].dot(faces.normal[iface3]) - 1.0);
            double d3 = std::abs(faces.normal[iface1].dot(faces.normal[iface4]) - 1.0);
            if (d1 + d2 + d3 > 1.0e-5) {
                std::cout << "\tSubfaces are not co-directed (" + side_to_string(s) + " side)\n";
                print_info(ic);
                return -1;
            }
        }
    }

    return 0;
}

int SoaCell::check_connectivity(int ic) const {
    if (is_undefined(ic)) {
        return 0;
    }

    // Через обычные грани существуют соседи
    for (int iface = face_begin[ic]; iface < face_begin[ic + 1]; ++iface) {
        if (faces.is_undefined(iface)) {
            // Проверим обнуление параметров
            if (faces.adjacent[iface].rank >= 0) {
                std::cout << "\tUndefined face, adjacent.rank >= 0\n";
                print_info(ic);
                return -1;
            }
            if (faces.adjacent[iface].alien >= 0) {
                std::cout << "\tUndefined face, adjacent.alien >= 0\n";
                print_info(ic);
                return -1;
            }
            if (faces.adjacent[iface].index >= 0) {
                std::cout << "\tUndefined face, adjacent.index >= 0\n";
                print_info(ic);
                return -1;
            }
            continue;
        }

        // Простая граничная грань
        if (faces.boundary[iface] != Boundary::ORDINARY &&
            faces.boundary[iface] != Boundary::PERIODIC) {

            // Грань должна ссылаться на искомую ячейку
            if (faces.adjacent[iface].rank != rank[ic] ||
                faces.adjacent[iface].index != index[ic] ||
                faces.adjacent[iface].alien >= 0) {
                std::cout << "\tBoundary face should point to origin cell\n";
                print_info(ic);
                return -1;
            }

            continue;
        }

        index_t jc = -1;
        auto& adj = faces.adjacent[iface];

        if (adj.rank < 0) {
            std::cout << "\tadjacent.rank < 0\n";
            print_info(ic);
            return -1;
        }

        if (adj.rank == mpi::rank()) {
            // Локальная ячейка
            if (adj.alien >= 0) {
                std::cout << "\tadjacent.alien >= 0\n";
                print_info(ic);
                return -1;
            }
            if (adj.index < 0 || adj.index >= size()) {
                std::cout << "\tLocal neighbor out of range\n";
                print_info(ic);
                return -1;
            }
            if (adj.index == ic) {
                std::cout << "\tSelf reference\n";
                print_info(ic);
                return -1;
            }
            jc = adj.index;
        }
        else {
            // Удаленная ячейка
            if (adj.index < 0) {
                std::cout << "\tRemote neighbor wrong index\n";
                print_info(ic);
                return -1;
            }
            /*  HOW TO HANDLE THIS????
            if (adj.alien < 0 || adj.alien >= aliens.size()) {
                std::cout << "\tRemote neighbor out of range\n";
                cell.print_info();
                return -1;
            }
            neib_ptr = &aliens[adj.alien];
            */
        }

        if (jc < 0) {
            std::cout << "\tUndefined neighbor\n";
            print_info(ic);
            return -1;
        }

        Vector3d fc(0.0, 0.0, 0.0);
        for (int i = 0; i < VpF(dim); ++i) {
            fc += verts[node_begin[ic] + faces.vertices[iface][i]];
        }
        fc /= VpF(dim);

        // Сосед должен иметь точно такую же грань, но с противоположной нормалью,
        // и ссылаться на текущую ячейку
        int counter = 0;

        for (index_t jface = face_begin[jc]; jface < face_begin[jc + 1]; ++jface) {
            if (faces.is_undefined(jface)) {
                continue;
            }

            if (faces.boundary[jface] != Boundary::ORDINARY &&
                faces.boundary[jface] != Boundary::PERIODIC) {
                // простая граничная грань
                continue;
            }

            Vector3d nfc(0.0, 0.0, 0.0);
            for (int i = 0; i < VpF(dim); ++i) {
                nfc += verts[node_begin[jc] + faces.vertices[jface][i]];
            }
            nfc /= VpF(dim);

            if ((fc - nfc).norm() > 1.0e-6 * linear_size(ic)) {
                continue;
            }

            // мы нашли соответствующую грань соседа
            ++counter;

            // нормали противоположны
            if (std::abs(faces.normal[iface].dot(faces.normal[jface]) + 1.0) > 1.0e-6) {
                std::cout << "\tOpposite faces have not opposite normals\n";
                std::cout << "\tCurrent cell:\n";
                print_info(ic);
                std::cout << "\tNeighbor:\n";
                print_info(jc);
                return -1;
            }
            // площади совпадают
            if (std::abs(faces.area[ic] - faces.area[jface]) > 1.0e-6 * linear_size(ic)) {
                std::cout << "\tOpposite faces have different area\n";
                std::cout << "\tCurrent cell:\n";
                print_info(ic);
                std::cout << "\tNeighbor:\n";
                print_info(jc);
                return -1;
            }
            // Указывает на исходную ячейку
            if (faces.adjacent[jface].index != ic) {
                std::cout << "\tWrong connection (index != ic)\n";
                std::cout << "\tCurrent cell:\n";
                print_info(ic);
                std::cout << "\tNeighbor:\n";
                print_info(jc);
                return -1;
            }
            // Ранг смежной (исходная) больше нуля
            if (faces.adjacent[jface].rank < 0) {
                std::cout << "\tWrong adjacent (rank < 0)\n";
                std::cout << "\tCurrent cell:\n";
                print_info(ic);
                std::cout << "\tNeighbor:\n";
                print_info(jc);
                return -1;
            }
            // Ранг смежной (исходная) равен рангу процесса
            if (faces.adjacent[jface].rank != mpi::rank()) {
                std::cout << "\tWrong adjacent (rank)\n";
                std::cout << "\tCurrent cell:\n";
                print_info(ic);
                std::cout << "\tNeighbor:\n";
                print_info(jc);
                return -1;
            }
            if (adj.rank == mpi::rank()) {
                // Локальный сосед
                if (faces.adjacent[jface].alien >= 0) {
                    std::cout << "\tWrong connection (alien >= 0)\n";
                    std::cout << "\tCurrent cell:\n";
                    print_info(ic);
                    std::cout << "\tNeighbor:\n";
                    print_info(jc);
                    return -1;
                }
            }
            else {
                // Удаленный сосед
                if (faces.adjacent[jface].alien < 0) {
                    // Проверка не актуальна сразу после построения списка соседей,
                    // face.alien для alien ячеек не получены, да и не важны.
                    /*
                    std::cout << "\tWrong connection\n";
                    std::cout << "\tCurrent cell:\n";
                    cell.print_info();
                    std::cout << "\tNeighbor:\n";
                    neib.print_info();
                    return -1;
                     */
                }
            }
        }
        if (counter < 1) {
            std::cout << "\tHas no neighbor across ordinary " << side_to_string(iface) << ")\n";
            print_info(ic);
            std::cout << "\tNeighbor:\n";
            print_info(faces.adjacent[iface].index);
            return -1;
        }
        if (counter > 1) {
            std::cout << "\tMore than one neighbor across ordinary face\n";
            print_info(ic);
            return -1;
        }
    }

    return 0;
}

int SoaMesh::check_base() {
    if (cells.empty()) {
        if (mpi::single()) {
            std::cout << "\tEmpty storage\n";
            return -1;
        } else {
            return 0;
        }
    }

    auto dim = cells.dim;

    if (dim != 2 && dim != 3) {
        std::cout << "\tDimension is not 2 or 3\n";
        return -1;
    }

    int res = 0;
    for (int ic = 0; ic < cells.size(); ++ic) {
        if (cells.index[ic] < 0 || cells.index[ic] != ic) {
            std::cout << "\tWrong cell index\n";
            return -1;
        }

        if (cells.rank[ic] < 0 || cells.rank[ic] != mpi::rank()) {
            std::cout << "\tWrong cell rank\n";
            return -1;
        }

        // Число граней
        for (int i = 0; i < FpC(dim); ++i) {
            if (cells.faces.is_undefined(cells.face_begin[ic] + i)) {
                std::cout << "\tCell has no one of main faces\n";
                cells.print_info(ic);
                return -1;
            }
        }
        if (cells.face_count(ic) > FpC(dim)) {
            std::cout << "\tCell has too much faces (" << cells.face_count(ic) << ")\n";
            cells.print_info(ic);
            return -1;
        }

        // Число вершин ???

        // Правильное задание геометрии
        res = cells.check_geometry(ic);
        if (res < 0) return res;

        // Грани правльно ориентированы
        res = cells.check_base_face_orientation(ic);
        if (res < 0) return res;

        // Порядок основных вершин
        res = cells.check_base_vertices_order(ic);
        if (res < 0) return res;

        // Проверка смежности
        res = cells.check_connectivity(ic);
        if (res < 0) return res;
    }

    return 0;
}

int SoaMesh::check_refined() {
    if (cells.empty()) {
        if (mpi::single()) {
            std::cout << "\tEmpty storage\n";
            return -1;
        } else {
            return 0;
        }
    }

    auto dim = cells.dim;

    if (dim != 2 && dim != 3) {
        std::cout << "\tDimension is not 2 or 3\n";
        return -1;
    }

    int res = 0;
    for (int ic = 0; ic < cells.size(); ++ic) {
        if (cells.is_undefined(ic)) {
            std::cout << "\tUndefined cell\n";
            return -1;
        }

        if (cells.index[ic] < 0 || cells.index[ic] != ic) {
            std::cout << "\tWrong cell index\n";
            return -1;
        }

        if (cells.rank[ic] < 0 || cells.rank[ic] != mpi::rank()) {
            std::cout << "\tWrong cell rank\n";
            return -1;
        }

        // Число граней
        for (int i = 0; i < FpC(dim); ++i) {
            if (cells.faces.is_undefined(cells.face_begin[ic] + i)) {
                std::cout << "\tCell has no one of main faces\n";
                cells.print_info(ic);
                return -1;
            }
        }

        // Вершины дублируются
        for (int i = cells.node_begin[ic]; i < cells.node_begin[ic + 1]; ++i) {
            for (int j = i + 1; j < cells.node_begin[ic + 1]; ++j) {
                double dist = (cells.verts[i] - cells.verts[j]).norm();
                if (dist < 1.0e-5 * cells.linear_size(ic)) {
                    std::cout << "\tIdentical vertices\n";
                    cells.print_info(ic);
                    return -1;
                }
            }
        }

        // Правильное задание геометрии
        res = cells.check_geometry(ic);
        if (res < 0) return res;

        // Грани правльно ориентированы
        res = cells.check_base_face_orientation(ic);
        if (res < 0) return res;

        // Порядок основных вершин
        res = cells.check_base_vertices_order(ic);
        if (res < 0) return res;

        // Проверка сложных граней
        res = cells.check_complex_faces(ic);
        if (res < 0) return res;

        // Проверка смежности
        res = cells.check_connectivity(ic);
        if (res < 0) return res;
    }

    return 0;
}

} // namespace zephyr::mesh