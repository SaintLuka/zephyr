#include <iostream>
#include <algorithm>
#include <format.h>

#include <zephyr/geom/boundary.h>
#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/side.h>
#include <zephyr/geom/indexing.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/primitives/quad.h>
#include <zephyr/utils/json.h>
#include <zephyr/mesh/raw/raw_cells.h>

namespace zephyr::geom::generator {

using namespace mesh;

Rectangle::Rectangle(const Json& config)
    : Generator("rectangle"),
      x_min_(0.0), x_max_(1.0),
      y_min_(0.0), y_max_(1.0) {

    axial_ = false;
    if (config["axial"]) {
        axial_ = config["axial"].as<bool>();
    }

    if (!config["geometry"]) {
        throw std::runtime_error("Rectangle config doesn't contain key 'geometry'");
    }

    x_min_ = config["geometry"]["x_min"].as<double>();
    x_max_ = config["geometry"]["x_max"].as<double>();
    y_min_ = config["geometry"]["y_min"].as<double>();
    y_max_ = config["geometry"]["y_max"].as<double>();

    if (!config["bounds"]) {
        throw std::runtime_error("Rectangle config doesn't contain key 'bounds'");
    }
    bounds_.left   = boundary_from_string(config["bounds"]["left"].as<std::string>());
    bounds_.right  = boundary_from_string(config["bounds"]["right"].as<std::string>());
    bounds_.bottom = boundary_from_string(config["bounds"]["bottom"].as<std::string>());
    bounds_.top    = boundary_from_string(config["bounds"]["top"].as<std::string>());

    if (config["voronoi"]) {
        voronoi_ = config["voronoi"].as<bool>();
    }

    // Адаптивная по умолчанию
    adaptive_ = !voronoi_;

    if (!config["size"]) {
        throw std::runtime_error("Rectangle config doesn't contain key 'size'");
    }

    if (config["size"].is_number()) {
        set_size(config["cells"].as<int>());
    } else {
        int nx(0), ny(0);
        if (config["size"]["nx"]) {
            nx = config["size"]["nx"].as<int>();
        }
        if (config["size"]["ny"]) {
            ny = config["size"]["ny"].as<int>();
        }
        if (nx > 0 && ny > 0) {
            set_sizes(nx, ny);
        } else if (nx > 0) {
            set_nx(nx);
        } else if (ny > 0) {
            set_ny(ny);
        } else {
            throw std::runtime_error(std::format("Rectangle config strange sizes: {}, {}, "
                                                 "setup key size.nx, size.ny or both", nx, ny));
        }
    }
}

Rectangle::Rectangle()
    : Rectangle(0.0, 1.0, 0.0, 1.0, false) {
    set_boundaries(Boundaries{});
}

Rectangle::Rectangle(double x_min, double x_max, double y_min, double y_max, bool voronoi) :
        Generator("rectangle"),
        x_min_(x_min), x_max_(x_max),
        y_min_(y_min), y_max_(y_max),
        voronoi_(voronoi) {
    if (!voronoi) {
        set_adaptive(true);
    }
    check_params();
}

Box Rectangle::bbox() const {
    Vector3d v_min(x_min_, y_min_, 0.0);
    Vector3d v_max(x_max_, y_max_, 0.0);

    return {v_min, v_max};
}

void Rectangle::set_axial(bool axial) {
    axial_ = axial;

    if (axial_ && y_min_ == 0.0) {
        bounds_.bottom = Boundary::WALL;
    }
}

void Rectangle::set_nx(int nx) {
    if (nx < 1) {
        throw std::runtime_error("Rectangle: Nx < 1");
    }
    if (!voronoi_) {
        nx_ = nx;
        ny_ = std::max(int(std::round(nx_ * (y_max_ - y_min_) / (x_max_ - x_min_))), 1);
    }
    else {
        nx_ = std::max(int(std::round(nx * std::pow(3.0, 0.25)) / std::sqrt(2.0)), 1);
        ny_ = std::max(int(std::round(nx_ * (y_max_ - y_min_) / (x_max_ - x_min_) / std::sqrt(3.0))), 1);
    }
    compute_size();
}

void Rectangle::set_ny(int ny) {
    if (ny < 1) {
        throw std::runtime_error("Rectangle: Ny < 1");
    }
    if (!voronoi_) {
        ny_ = ny;
        nx_ = std::max(int(std::round(ny * (x_max_ - x_min_) / (y_max_ - y_min_))), 1);
    }
    else {
        ny_ = std::max(int(std::round(ny / std::pow(3.0, 0.25)) / std::sqrt(2.0)), 1);
        nx_ = std::max(int(std::round(ny_ * (x_max_ - x_min_) / (y_max_ - y_min_) * std::sqrt(3.0))), 1);
    }
    compute_size();
}

void Rectangle::set_sizes(int nx, int ny) {
    if (nx < 1 || ny < 1) {
        throw std::runtime_error("Rectangle: Nx < 1 or Ny < 1");
    }
    if (!voronoi_) {
        nx_ = nx;
        ny_ = ny;
    }
    else {
        nx_ = std::max(int(std::round(nx * std::pow(3.0, 0.25)) / std::sqrt(2.0)), 1);
        ny_ = std::max(int(std::round(ny / std::pow(3.0, 0.25)) / std::sqrt(2.0)), 1);
    }
    compute_size();

    double dx = (x_max_ - x_min_) / nx_;
    double dy = (y_max_ - y_min_) / ny_;

    double dmax = std::max(dx, dy);
    double dmin = std::min(dx, dy);
    if (dmax / dmin > 2.0) {
        std::cerr << "Rectangle Warning: Large aspect ratio (> 2)\n";
    }
    if (dmax / dmin > 1.0e3) {
        std::cerr << "Rectangle Warning: Huge aspect ratio (> 1000)\n";
    }
}

void Rectangle::set_size(int N) {
    if (N < 1) {
        std::cerr << "Rectangle Error: N < 1\n";
        throw std::runtime_error("Rectangle Error: N < 1");
    }

    if (!voronoi_) {
        double d = std::sqrt((x_max_ - x_min_) * (y_max_ - y_min_) / N);
        nx_ = std::max(int(std::round((x_max_ - x_min_) / d)), 1);
        ny_ = std::max(int(std::round((y_max_ - y_min_) / d)), 1);
    }
    else {
        double a = std::sqrt(2 * (x_max_ - x_min_) * (y_max_ - y_min_) / (std::sqrt(3.0) * N));
        double h = std::sqrt(3.0) * a;

        nx_ = std::max(int(std::round((x_max_ - x_min_) / a)), 1);
        ny_ = std::max(int(std::round((y_max_ - y_min_) / h)), 1);
    }
    compute_size();
}

void Rectangle::set_boundaries(Boundaries bounds) {
    bounds_ = bounds;
    if (periodic_along_x()) {
        bounds_.left = bounds_.right = Boundary::PERIODIC;
    }
    if (periodic_along_y()) {
        bounds_.bottom = bounds_.top = Boundary::PERIODIC;
    }
    if (axial_ && y_min_ == 0.0) {
        bounds_.bottom = Boundary::WALL;
    }
}

void Rectangle::set_adaptive(bool adaptive) {
    if (!voronoi_) {
        adaptive_ = adaptive;
    }
    else {
        adaptive_ = false;
    }
}

double Rectangle::x_min() const {
    return x_min_;
}

double Rectangle::x_max() const {
    return x_max_;
}

double Rectangle::y_min() const {
    return y_min_;
}

double Rectangle::y_max() const {
    return y_max_;
}

int Rectangle::nx() const {
    return nx_;
}

int Rectangle::ny() const {
    return ny_;
}

Rectangle::Boundaries Rectangle::bounds() const {
    return bounds_;
}

bool Rectangle::periodic_along_x() const {
    return bounds_.left == Boundary::PERIODIC || bounds_.right == Boundary::PERIODIC;
}

bool Rectangle::periodic_along_y() const {
    return bounds_.bottom == Boundary::PERIODIC || bounds_.top == Boundary::PERIODIC;
}

void Rectangle::check_params() const {
    if (x_min_ >= x_max_) {
        throw std::runtime_error("Rectangle::check_params: x_min >= x_max");
    }
    if (y_min_ >= y_max_) {
        throw std::runtime_error("Rectangle::check_params: y_min >= y_max");
    }
}

void Rectangle::compute_size() {
    if (!voronoi_) {
        size_ = nx_ * ny_;
    } else {
        size_ = 2 * nx_ * ny_;
    }
    if (size_ > max_grid_size) {
        throw std::runtime_error(std::format("Generator::check_size: attempt to create mesh "
                                             "that contains more than {} elements", max_grid_size));
    }
}

Grid Rectangle::make() const {
    if (voronoi_) {
        return create_voronoi();
    }
    if (adaptive_) {
        return create_classic_amr();
    }
    return create_classic();
}

Grid Rectangle::create_classic() const {
    check_size(size_);

    double dx = (x_max_ - x_min_) / nx_;
    double dy = (y_max_ - y_min_) / ny_;

    Grid grid;
    grid.reserve_nodes((nx_ + 1) * (ny_ + 1));
    grid.reserve_cells(nx_ * ny_);

    std::vector nodes(nx_ + 1, std::vector<GNode::Ptr>(ny_ + 1));
    for (int i = 0; i <= nx_; ++i) {
        for (int j = 0; j <= ny_; ++j) {
            double x = x_min_ + i * dx;
            double y = y_min_ + j * dy;
            nodes[i][j] = GNode::create({x, y, 0.0});
        }
    }

    // Ячейки как полигоны, обход граней и вершин против часовой
    // стрелки, начиная с нижней левой вершины (нижней грани)
    std::vector<Boundary> bc(4);
    std::vector<GNode::Ptr> quad_nodes(4);
    for (int i = 0; i < nx_; ++i) {
        bc[Side2D::L] = i == 0      ? bounds_.left   : Boundary::INNER;
        bc[Side2D::R] = i == nx_-1 ? bounds_.right  : Boundary::INNER;
        for (int j = 0; j < ny_; ++j) {
            bc[Side2D::B] = j == 0      ? bounds_.bottom : Boundary::INNER;
            bc[Side2D::T] = j == ny_-1 ? bounds_.top    : Boundary::INNER;

            using indexing::quad::vs;
            quad_nodes[vs<0, 0>()] = nodes[i][j];
            quad_nodes[vs<0, 1>()] = nodes[i][j+1];
            quad_nodes[vs<1, 0>()] = nodes[i+1][j];
            quad_nodes[vs<1, 1>()] = nodes[i+1][j+1];

            grid.add_cell(CellType::QUAD, quad_nodes, bc);
        }
    }
    return grid;
}

Grid Rectangle::create_classic_amr() const {
    check_size(size_);

    double dx = (x_max_ - x_min_) / nx_;
    double dy = (y_max_ - y_min_) / ny_;

    Grid grid;
    grid.reserve_nodes((2 * nx_ + 1) * (2 * ny_ + 1));
    grid.reserve_cells(nx_ * ny_);

    std::vector nodes(2 * nx_ + 1, std::vector<GNode::Ptr>(2 * ny_ + 1));
    for (int i = 0; i <= 2 * nx_; ++i) {
        for (int j = 0; j <= 2 * ny_; ++j) {
            double x = x_min_ + 0.5 * i * dx;
            double y = y_min_ + 0.5 * j * dy;
            nodes[i][j] = GNode::create({x, y, 0.0});
        }
    }

    // AMR-ячейки, обход граней как в Side2D, вершины в Z-порядке
    std::vector<Boundary> bc(4);
    for (int i = 0; i < nx_; ++i) {
        bc[Side2D::L] = i == 0      ? bounds_.left   : Boundary::INNER;
        bc[Side2D::R] = i == nx_-1 ? bounds_.right  : Boundary::INNER;
        for (int j = 0; j < ny_; ++j) {
            bc[Side2D::B] = j == 0      ? bounds_.bottom : Boundary::INNER;
            bc[Side2D::T] = j == ny_-1 ? bounds_.top    : Boundary::INNER;
            grid.add_cell(
                CellType::AMR2D, {
                    nodes[2*i][2*j + 0], nodes[2*i + 1][2*j + 0], nodes[2*i + 2][2*j + 0],
                    nodes[2*i][2*j + 1], nodes[2*i + 1][2*j + 1], nodes[2*i + 2][2*j + 1],
                    nodes[2*i][2*j + 2], nodes[2*i + 1][2*j + 2], nodes[2*i + 2][2*j + 2],
                }, bc);
        }
    }
    return grid;
}

Grid Rectangle::create_voronoi() const {
    check_size(size_);

    size_t Nb = nx_ * ny_;
    double DX = x_max_ - x_min_;
    double DY = y_max_ - y_min_;
    size_t Ny = size_t(std::floor(std::sqrt(std::sqrt(3.0) * DY * Nb / 2.0 / DX))) + 1;
    double h = DY / 2.0 / Ny;
    double D = h / std::sqrt(0.75);

    size_t Nx = size_t(std::floor(DX / 1.5 / D));

    double x_shift = x_min_ + (DX - Nx * 1.5 * D) / 2.0 - 0.5 * D;
    double y_shift = y_min_;

    // Вершины в виде таблицы
    std::vector vertices(Nx + 2, std::vector<GNode::Ptr>(2 * Ny + 1, nullptr));

    for (size_t j = 0; j <= 2 * Ny; ++j) {
        double y = y_shift + h * j;

        // Часть вершин на левой границе пропускаем
        if (j % 2 == 0) {
            vertices[0][j] = GNode::create({x_min_, y, 0.0});
            vertices[0][j]->bc = bounds_.left;
        }

        for (size_t i = 1; i <= Nx; ++i) {
            double x = x_shift + double(3 * i - (i + j) % 2) * 0.5 * D;
            vertices[i][j] = GNode::create({x, y, 0.0});

            if (j == 0) {
                vertices[i][j]->bc = bounds_.bottom;
            } else if (j == 2 * Ny) {
                vertices[i][j]->bc = bounds_.top;
            }
        }

        // Часть вершин на правой границе пропускаем
        if (j == 0 || j == 2 * Ny || j % 2 == Nx % 2) {
            vertices[Nx + 1][j] = GNode::create({x_max_, y, 0.0});
            vertices[Nx + 1][j]->bc = bounds_.right;
        }
    }

    using VList = std::vector<GNode::Ptr>;

    auto erase_nans = [](VList& vlist) {
        const auto to_remove = std::ranges::remove_if(vlist,
            [](GNode::Ref v) -> bool { return !v; }).begin();
        vlist.erase(to_remove, vlist.end());
    };

    const double eps = 1.0e-9;

    auto get_bounds = [this, eps](VList& vlist) -> std::vector<Boundary> {
        std::vector bounds(vlist.size(), Boundary::INNER);
        for (int i = 0; i < vlist.size(); ++i) {
            int j = (i + 1) % vlist.size();
            Vector3d v1 = vlist[i]->pos;
            Vector3d v2 = vlist[j]->pos;
            if (std::max(v1.x(), v2.x()) <= x_min_ + eps) { bounds[i] = bounds_.left; }
            if (std::min(v1.x(), v2.x()) >= x_max_ - eps) { bounds[i] = bounds_.right; }
            if (std::max(v1.y(), v2.y()) <= y_min_ + eps) { bounds[i] = bounds_.bottom; }
            if (std::min(v1.y(), v2.y()) >= y_max_ - eps) { bounds[i] = bounds_.top; }
        }
        return bounds;
    };

    Grid grid;
    grid.reserve_nodes((Nx + 2)*(2 * Ny + 1));
    grid.reserve_cells(size_);
    for (size_t i = 0; i <= Nx; ++i) {
        if (i % 2 == 0) {
            for (size_t j = 0; j < 2 * Ny - 1; j += 2) {
                VList vlist = {
                    vertices[i][j], vertices[i + 1][j],
                    vertices[i + 1][j + 1], vertices[i + 1][j + 2],
                    vertices[i][j + 2], vertices[i][j + 1]
                };
                erase_nans(vlist);
                grid.add_cell(CellType::POLYGON, vlist, get_bounds(vlist));
            }
        } else {
            // j == 0
            VList vlist = {
                vertices[i][0], vertices[i + 1][0],
                vertices[i + 1][1], vertices[i][1]
            };
            grid.add_cell(CellType::POLYGON, vlist, get_bounds(vlist));

            for (size_t j = 1; j < 2 * Ny - 2; j += 2) {
                vlist = {
                    vertices[i][j], vertices[i + 1][j],
                    vertices[i + 1][j + 1], vertices[i + 1][j + 2],
                    vertices[i][j + 2], vertices[i][j + 1]
                };
                erase_nans(vlist);
                grid.add_cell(CellType::POLYGON, vlist, get_bounds(vlist));
            }

            // j == last; i % 2 == 1
            vlist = {
                vertices[i][2 * Ny - 1], vertices[i + 1][2 * Ny - 1],
                vertices[i + 1][2 * Ny], vertices[i][2 * Ny]
            };
            grid.add_cell(CellType::POLYGON, vlist, get_bounds(vlist));
        }
    }

    return grid;
}

RawCells Rectangle::make_cells(bool unique_nodes) const {
    if (voronoi_ && !adaptive_) {
        throw std::runtime_error("Rectangle::initialize: can't initialize voronoi grid and classic cartesian");
    }
    check_size(size_);

    bool x_period = periodic_along_x();
    bool y_period = periodic_along_y();

    double hx = (x_max_ - x_min_) / nx_;
    double hy = (y_max_ - y_min_) / ny_;

    auto get_index = [this](index_t i, index_t j) -> index_t {
        return i * ny_ + j;
    };

    auto get_index_pair = [this](index_t n) -> std::array<index_t, 2> {
        return {n / ny_, n % ny_};
    };

    auto get_vertex = [=, this](index_t i, index_t j) -> Vector3d {
        return {
                x_min_ + ((x_max_ - x_min_) * i) / nx_,
                y_min_ + ((y_max_ - y_min_) * j) / ny_,
                0.0
        };
    };

    auto neib_index = [=, this](index_t i, index_t j, Side2D side) -> index_t {
        if (side == Side2D::LEFT) {
            return i == 0 && !x_period ?  get_index(i, j) : get_index((i - 1 + nx_) % nx_, j);
        }
        if (side == Side2D::RIGHT) {
            return i == nx_ - 1 && !x_period ? get_index(i, j) : get_index((i + 1 + nx_) % nx_, j);
        }
        if (side == Side2D::BOTTOM) {
            return j == 0 && !y_period ? get_index(i, j) : get_index(i, (j - 1 + ny_) % ny_);
        }
        if (side == Side2D::TOP) {
            return j == ny_ - 1 && !y_period ? get_index(i, j): get_index(i, (j + 1 + ny_) % ny_);
        }
        throw std::runtime_error("Strange side #142");
    };

    RawCells cells({
        .dim = 2,
        .adaptive = true,
        .linear = true,
        .axial = axial_,
        .nodes = unique_nodes
    });

    cells.resize_amr(size_);

    static constexpr int n_faces = 8;
    static constexpr int n_nodes = 9;

    for (index_t ic = 0; ic < size_; ++ic) {
        auto[i, j] = get_index_pair(ic);

        cells.next[ic] = ic;
        cells.rank[ic] = 0;
        cells.index[ic] = ic;

        cells.flag[ic] = 0;
        cells.level[ic] = 0;
        cells.b_idx[ic] = ic;
        cells.z_idx[ic] = 0;

        SqQuad quad(
            get_vertex(i, j),
            get_vertex(i + 1, j),
            get_vertex(i, j + 1),
            get_vertex(i + 1, j + 1));

        cells.center[ic] = quad.vs<0, 0>();
        cells.volume[ic] = hx * hy;
        cells.volume_alt[ic] = NAN;
        cells.verts.offsets[ic] = ic * n_nodes;
        cells.verts.offsets[ic + 1] = (ic + 1) * n_nodes;

        // INIT FACES
        cells.faces.offsets[ic] = ic * n_faces;
        cells.faces.offsets[ic + 1] = (ic + 1) * n_faces;

        for (auto iface: cells.faces.range(ic)) {
            cells.faces.set_undefined(iface);
            cells.faces.area_alt[iface] = NAN;
        }

        index_t iface = ic * n_faces;

        cells.faces.boundary[iface + Side2D::L] = i > 0 ? Boundary::INNER : bounds_.left;
        cells.faces.boundary[iface + Side2D::R] = i < nx_ - 1 ? Boundary::INNER : bounds_.right;
        cells.faces.boundary[iface + Side2D::B] = j > 0 ? Boundary::INNER : bounds_.bottom;
        cells.faces.boundary[iface + Side2D::T] = j < ny_ - 1 ? Boundary::INNER : bounds_.top;

        for (auto side: Side2D::items()) {
            cells.faces.adjacent.rank[iface + side] = 0;
            cells.faces.adjacent.index[iface + side] = neib_index(i, j, side);
            cells.faces.adjacent.ghost[iface + side] = -1;
            cells.faces.adjacent.basic[iface + side] = ic;
            cells.faces.vertices[iface + side].fill(-1);
        }

        cells.faces.normal[iface + Side2D::L] = -Vector3d::UnitX();
        cells.faces.normal[iface + Side2D::R] =  Vector3d::UnitX();
        cells.faces.normal[iface + Side2D::B] = -Vector3d::UnitY();
        cells.faces.normal[iface + Side2D::T] =  Vector3d::UnitY();

        cells.faces.center[iface + Side2D::L] = quad.vs<-1, 0>();
        cells.faces.center[iface + Side2D::R] = quad.vs<+1, 0>();
        cells.faces.center[iface + Side2D::B] = quad.vs<0, -1>();
        cells.faces.center[iface + Side2D::T] = quad.vs<0, +1>();

        cells.faces.area[iface + Side2D::L] = hy;
        cells.faces.area[iface + Side2D::R] = hy;
        cells.faces.area[iface + Side2D::B] = hx;
        cells.faces.area[iface + Side2D::T] = hx;

        cells.faces.vertices[iface + Side2D::L] = indexing::amr::sf(Side2D::L);
        cells.faces.vertices[iface + Side2D::R] = indexing::amr::sf(Side2D::R);
        cells.faces.vertices[iface + Side2D::B] = indexing::amr::sf(Side2D::B);
        cells.faces.vertices[iface + Side2D::T] = indexing::amr::sf(Side2D::T);

        for (index_t jn = 0; jn < n_nodes; ++jn) {
            cells.verts[ic * n_nodes + jn] = quad[jn];
        }

        if (axial_) {
            // "Альтернативный" объем ячейки и площади граней
            cells.volume_alt[ic] = hx * hy * quad.vs<0, 0>().y();
            cells.faces.area_alt[iface + Side2D::L] = hy * quad.vs<-1, 0>().y();
            cells.faces.area_alt[iface + Side2D::R] = hy * quad.vs<+1, 0>().y();
            cells.faces.area_alt[iface + Side2D::B] = hx * quad.vs< 0,-1>().y();
            cells.faces.area_alt[iface + Side2D::T] = hx * quad.vs< 0,+1>().y();

            // Смещения барицентров, есть необходимость?
            cells.center[ic].y() += hy*hy / (12.0 * quad.vs<0, 0>().y());
            cells.faces.center[iface + Side2D::L].y() += hy*hy / (12.0 * quad.vs<-1, 0>().y());
            cells.faces.center[iface + Side2D::R].y() += hy*hy / (12.0 * quad.vs<+1, 0>().y());
        }
    }
    return cells;
}

} // namespace zephyr::geom::generator