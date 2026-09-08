#include <iostream>
#include <format.h>
#include <algorithm>

#include <zephyr/geom/side.h>
#include <zephyr/geom/indexing.h>
#include <zephyr/geom/boundary.h>
#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/primitives/cube.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/utils/json.h>
#include <zephyr/mesh/euler/amr_cells.h>

namespace zephyr::geom::generator {

using namespace mesh;

Cuboid::Cuboid(const Json& config)
    : Generator("cuboid"),
      x_min_(0.0), x_max_(1.0),
      y_min_(0.0), y_max_(1.0),
      z_min_(0.0), z_max_(1.0) {

    // Адаптивная по умолчанию
    adaptive_ = true;

    if (!config["geometry"]) {
        throw std::runtime_error("Cuboid config doesn't contain key 'geometry'");
    }

    x_min_ = config["geometry"]["x_min"].as<double>();
    x_max_ = config["geometry"]["x_max"].as<double>();
    y_min_ = config["geometry"]["y_min"].as<double>();
    y_max_ = config["geometry"]["y_max"].as<double>();
    z_min_ = config["geometry"]["z_min"].as<double>();
    z_max_ = config["geometry"]["z_max"].as<double>();

    if (!config["bounds"]) {
        throw std::runtime_error("Cuboid config doesn't contain key 'bounds'");
    }
    bounds_.left   = boundary_from_string(config["bounds"]["left"].as<std::string>());
    bounds_.right  = boundary_from_string(config["bounds"]["right"].as<std::string>());
    bounds_.bottom = boundary_from_string(config["bounds"]["bottom"].as<std::string>());
    bounds_.top    = boundary_from_string(config["bounds"]["top"].as<std::string>());
    bounds_.back   = boundary_from_string(config["bounds"]["back"].as<std::string>());
    bounds_.front  = boundary_from_string(config["bounds"]["front"].as<std::string>());

    if (!config["size"]) {
        throw std::runtime_error("Cuboid config doesn't contain key 'size'");
    }

    if (config["size"].is_number()) {
        set_size(config["size"].as<int>());
    }
    else {
        int nx(0), ny(0), nz(0);
        if (config["size"]["nx"]) {
            nx = config["size"]["nx"].as<int>();
        }
        if (config["size"]["ny"]) {
            ny = config["size"]["ny"].as<int>();
        }
        if (config["size"]["nz"]) {
            nz = config["size"]["nz"].as<int>();
        }
        if (nx > 0 && ny > 0 && nz > 0) {
            set_sizes(nx, ny, nz);
        }
        else {
            if (nx > 0) {
                set_nx(nx);
            } else if (ny > 0) {
                set_ny(ny);
            } else if (nz > 0) {
                set_nz(nz);
            } else {
                throw std::runtime_error(std::format("Cuboid config strange sizes: {}, {}, {}, "
                                                     "setup key size.nx, size.ny, size.nz or all of them", nx, ny, nz));
            }
        }
    }
}

Cuboid::Cuboid(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max) :
        Generator("cuboid"),
        x_min_(x_min), x_max_(x_max),
        y_min_(y_min), y_max_(y_max),
        z_min_(z_min), z_max_(z_max) {
    // Адаптивная по умолчанию
    adaptive_ = true;
    check_params();
}

Box Cuboid::bbox() const {
    Vector3d vmin(x_min_, y_min_, z_min_);
    Vector3d vmax(x_max_, y_max_, z_max_);

    return {vmin, vmax};
}

void Cuboid::set_nx(int nx) {
    if (nx < 1) {
        throw std::runtime_error("Cuboid::set_nx: Nx < 1");
    }
    nx_ = nx;
    ny_ = int(round(nx * (y_max_ - y_min_) / (x_max_ - x_min_)));
    nz_ = int(round(nx * (z_max_ - z_min_) / (x_max_ - x_min_)));
    compute_size();
}

void Cuboid::set_ny(int ny) {
    if (ny < 1) {
        throw std::runtime_error("Cuboid::set_ny: Ny < 1");
    }
    ny_ = ny;
    nx_ = int(round(ny * (x_max_ - x_min_) / (y_max_ - y_min_)));
    nz_ = int(round(ny * (z_max_ - z_min_) / (y_max_ - y_min_)));
    compute_size();
}

void Cuboid::set_nz(int nz) {
    if (nz < 1) {
        throw std::runtime_error("Cuboid::set_nz: Nz < 1");
    }
    nz_ = nz;
    nx_ = int(round(nz * (x_max_ - x_min_) / (z_max_ - z_min_)));
    ny_ = int(round(nz * (y_max_ - y_min_) / (z_max_ - z_min_)));
    compute_size();
}

void Cuboid::set_sizes(int nx, int ny, int nz) {
    if (nx < 1 || ny < 1 || nz < 1) {
        throw std::runtime_error("Cuboid::set_sizes: Nx < 1 or Ny < 1 or Nz < 1");
    }
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
    compute_size();

    double dx = (x_max_ - x_min_) / nx_;
    double dy = (y_max_ - y_min_) / ny_;
    double dz = (z_max_ - z_min_) / nz_;

    double dmax = std::max(dx, std::max(dy, dz));
    double dmin = std::min(dx, std::min(dy, dz));
    if (dmax / dmin > 2.0) {
        std::cerr << "Cuboid Warning: Large aspect ratio (> 2)\n";
    }
    if (dmax / dmin > 1.0e3) {
        std::cerr << "Cuboid Warning: Huge aspect ratio (> 1000)\n";
    }
}

void Cuboid::set_size(int N) {
    double d = std::cbrt((x_max_ - x_min_) * (y_max_ - y_min_) * (z_max_ - z_min_) / N);
    nx_ = int(round((x_max_ - x_min_) / d));
    ny_ = int(round((y_max_ - y_min_) / d));
    nz_ = int(round((z_max_ - z_min_) / d));
    compute_size();
}

void Cuboid::set_boundaries(Boundaries bounds) {
    bounds_ = bounds;

    if (periodic_along_x()) {
        bounds_.left = bounds_.right = Boundary::PERIODIC;
    }
    if (periodic_along_y()) {
        bounds_.bottom = bounds_.top = Boundary::PERIODIC;
    }
    if (periodic_along_z()) {
        bounds_.back = bounds_.front = Boundary::PERIODIC;
    }
}

double Cuboid::x_min() const {
    return x_min_;
}

double Cuboid::x_max() const {
    return x_max_;
}

double Cuboid::y_min() const {
    return y_min_;
}

double Cuboid::y_max() const {
    return y_max_;
}

double Cuboid::z_min() const {
    return y_min_;
}

double Cuboid::z_max() const {
    return y_max_;
}

int Cuboid::nx() const {
    return nx_;
}

int Cuboid::ny() const {
    return ny_;
}

int Cuboid::nz() const {
    return ny_;
}

Cuboid::Boundaries Cuboid::bounds() const {
    return bounds_;
}

bool Cuboid::periodic_along_x() const {
    return bounds_.left == Boundary::PERIODIC || bounds_.right == Boundary::PERIODIC;
}

bool Cuboid::periodic_along_y() const {
    return bounds_.bottom == Boundary::PERIODIC || bounds_.top == Boundary::PERIODIC;
}

bool Cuboid::periodic_along_z() const {
    return bounds_.back == Boundary::PERIODIC || bounds_.front == Boundary::PERIODIC;
}

void Cuboid::check_params() const {
    if (x_min_ >= x_max_) {
        throw std::runtime_error("Cuboid::check_params: x_min >= x_max");
    }
    if (y_min_ >= y_max_) {
        throw std::runtime_error("Cuboid::check_params: y_min >= y_max");
    }
    if (z_min_ >= z_max_) {
        throw std::runtime_error("Cuboid::check_params: z_min >= z_max");
    }
    double dx = (x_max_ - x_min_) / nx_;
    double dy = (y_max_ - y_min_) / ny_;
    double dz = (y_max_ - y_min_) / ny_;

    double dmax = std::max(dx, std::max(dy, dz));
    double dmin = std::min(dx, std::min(dy, dz));
    if (dmax / dmin > 1.0e3) {
        std::cerr << "Cuboid::check_params warning: Huge aspect ratio (> 1000)\n";
    }
}

void Cuboid::compute_size() {
    size_ = nx_ * ny_ * nz_;
}

Grid Cuboid::make() const {
    check_size(size_);

    double dx = (x_max_ - x_min_) / nx_;
    double dy = (y_max_ - y_min_) / ny_;
    double dz = (z_max_ - z_min_) / nz_;

    Grid grid;
    grid.reserve_nodes((nx_ + 1) * (ny_ + 1) * (nz_ + 1));
    grid.reserve_cells(nx_ * ny_ * nz_);

    std::vector nodes(nx_ + 1, std::vector(ny_ + 1, std::vector<Node::Ptr>(nz_ + 1)));
    for (int i = 0; i <= nx_; ++i) {
        for (int j = 0; j <= ny_; ++j) {
            for (int k = 0; k <= nz_; ++k) {
                double x = x_min_ + i * dx;
                double y = y_min_ + j * dy;
                double z = z_min_ + k * dz;
                nodes[i][j][k] = Node::create({x, y, z});
            }
        }
    }

    std::vector<Boundary> bc(6);
    std::vector<Node::Ptr> cube_nodes(8);
    for (int i = 0; i < nx_; ++i) {
        bc[Side3D::L] = i == 0 ?      bounds_.left :  Boundary::INNER;
        bc[Side3D::R] = i == nx_-1 ? bounds_.right : Boundary::INNER;
        for (int j = 0; j < ny_; ++j) {
            bc[Side3D::B] = j == 0 ?      bounds_.bottom : Boundary::INNER;
            bc[Side3D::T] = j == ny_-1 ? bounds_.top : Boundary::INNER;
            for (int k = 0; k < nz_; ++k) {
                bc[Side3D::Z] = k == 0 ? bounds_.back : Boundary::INNER;
                bc[Side3D::F] = k == nz_-1 ? bounds_.front : Boundary::INNER;

                using indexing::hex::vs;
                cube_nodes[vs<0,0,0>()] = nodes[i][j][k];
                cube_nodes[vs<1,0,0>()] = nodes[i+1][j][k];
                cube_nodes[vs<1,1,0>()] = nodes[i+1][j+1][k];
                cube_nodes[vs<0,1,0>()] = nodes[i][j+1][k];
                cube_nodes[vs<0,0,1>()] = nodes[i][j][k+1];
                cube_nodes[vs<1,0,1>()] = nodes[i+1][j][k+1];
                cube_nodes[vs<1,1,1>()] = nodes[i+1][j+1][k+1];
                cube_nodes[vs<0,1,1>()] = nodes[i][j+1][k+1];
                grid.add_cell(CellType::HEXAHEDRON, cube_nodes, bc);
            }
        }
    }
    return grid;
}

AmrCells Cuboid::make_cells(bool unique_nodes) const {
    if (!adaptive_) {
        throw std::runtime_error("Cuboid::make_cells: can make_cells only AMR cartesian mesh");
    }

    bool x_period = periodic_along_x();
    bool y_period = periodic_along_y();
    bool z_period = periodic_along_z();

    double hx = (x_max_ - x_min_) / nx_;
    double hy = (y_max_ - y_min_) / ny_;
    double hz = (z_max_ - z_min_) / nz_;

    auto get_index = [=, this](index_t i, index_t j, index_t k) -> index_t {
        return nz_ * (ny_ * i + j) + k;
    };

    auto get_index_pair = [=, this](index_t n) -> std::array<index_t, 3> {
        return {(n / nz_) / ny_, (n / nz_) % ny_, n % nz_};
    };

    auto get_vertex = [=, this](index_t i, index_t j, index_t k) -> Vector3d {
        return {
                x_min_ + ((x_max_ - x_min_) * i) / nx_,
                y_min_ + ((y_max_ - y_min_) * j) / ny_,
                z_min_ + ((z_max_ - z_min_) * k) / nz_
        };
    };

    auto neib_index = [=, this](index_t i, index_t j, index_t k, Side3D side) -> index_t {
        if (side == Side3D::LEFT) {
            return i == 0 && !x_period ?  get_index(i, j, k) : get_index((i - 1 + nx_) % nx_, j, k);
        }
        if (side == Side3D::RIGHT) {
            return i == nx_ - 1 && !x_period ? get_index(i, j, k) : get_index((i + 1) % nx_, j, k);
        }
        if (side == Side3D::BOTTOM) {
            return j == 0 && !y_period ? get_index(i, j, k) : get_index(i, (j - 1 + ny_) % ny_, k);
        }
        if (side == Side3D::TOP) {
            return j == ny_ - 1 && !y_period ? get_index(i, j, k): get_index(i, (j + 1) % ny_, k);
        }
        if (side == Side3D::BACK) {
            return k == 0 && !z_period ? get_index(i, j, k) : get_index(i, j, (k - 1 + nz_) % nz_);
        }
        if (side == Side3D::FRONT) {
            return k == nz_ - 1 && !z_period ? get_index(i, j, k): get_index(i, j, (k + 1) % nz_);
        }
        throw std::runtime_error("Strange side #265");
    };

    AmrCells cells({
        .dim = 3,
        .adaptive = true,
        .linear = true,
        .axial = false,
        .nodes = unique_nodes
    });

    cells.resize_amr(size_);

    static constexpr int n_faces = 24;
    static constexpr int n_nodes = 27;

    for (index_t ic = 0; ic < size_; ++ic) {
        auto[i, j, k] = get_index_pair(ic);

        cells.next[ic] = ic;
        cells.rank[ic] = 0;
        cells.index[ic] = ic;

        cells.flag[ic] = 0;
        cells.level[ic] = 0;
        cells.b_idx[ic] = ic;
        cells.z_idx[ic] = 0;

        SqCube cube(get_vertex(i,   j,   k),
                    get_vertex(i+1, j,   k),
                    get_vertex(i,   j+1, k),
                    get_vertex(i+1, j+1, k),
                    get_vertex(i,   j,   k+1),
                    get_vertex(i+1, j,   k+1),
                    get_vertex(i,   j+1, k+1),
                    get_vertex(i+1, j+1, k+1));

        cells.center[ic] = cube.vs<0, 0, 0>();
        cells.volume[ic] = hx * hy * hz;
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

        cells.faces.boundary[iface + Side3D::L] = i > 0 ? Boundary::INNER : bounds_.left;
        cells.faces.boundary[iface + Side3D::R] = i < nx_ - 1 ? Boundary::INNER : bounds_.right;
        cells.faces.boundary[iface + Side3D::B] = j > 0 ? Boundary::INNER : bounds_.bottom;
        cells.faces.boundary[iface + Side3D::T] = j < ny_ - 1 ? Boundary::INNER : bounds_.top;
        cells.faces.boundary[iface + Side3D::Z] = k > 0 ? Boundary::INNER : bounds_.back;
        cells.faces.boundary[iface + Side3D::F] = k < nz_ - 1 ? Boundary::INNER : bounds_.front;

        for (auto side: Side3D::items()) {
            cells.faces.adjacent.rank[iface + side] = 0;
            cells.faces.adjacent.index[iface + side] = neib_index(i, j, k, side);
            cells.faces.adjacent.ghost[iface + side] = -1;
            cells.faces.adjacent.basic[iface + side] = ic;
            cells.faces.vertices[iface + side].fill(-1);
        }

        cells.faces.normal[iface + Side3D::L] = -Vector3d::UnitX();
        cells.faces.normal[iface + Side3D::R] =  Vector3d::UnitX();
        cells.faces.normal[iface + Side3D::B] = -Vector3d::UnitY();
        cells.faces.normal[iface + Side3D::T] =  Vector3d::UnitY();
        cells.faces.normal[iface + Side3D::Z] = -Vector3d::UnitZ();
        cells.faces.normal[iface + Side3D::F] =  Vector3d::UnitZ();

        cells.faces.center[iface + Side3D::L] = cube.vs<-1, 0, 0>();
        cells.faces.center[iface + Side3D::R] = cube.vs<+1, 0, 0>();
        cells.faces.center[iface + Side3D::B] = cube.vs< 0,-1, 0>();
        cells.faces.center[iface + Side3D::T] = cube.vs< 0,+1, 0>();
        cells.faces.center[iface + Side3D::Z] = cube.vs< 0, 0,-1>();
        cells.faces.center[iface + Side3D::F] = cube.vs< 0, 0,+1>();

        cells.faces.area[iface + Side3D::L] = hy * hz;
        cells.faces.area[iface + Side3D::R] = hy * hz;
        cells.faces.area[iface + Side3D::B] = hx * hz;
        cells.faces.area[iface + Side3D::T] = hx * hz;
        cells.faces.area[iface + Side3D::Z] = hx * hy;
        cells.faces.area[iface + Side3D::F] = hx * hy;

        cells.faces.vertices[iface + Side3D::L] = indexing::amr::sf(Side3D::L);
        cells.faces.vertices[iface + Side3D::R] = indexing::amr::sf(Side3D::R);
        cells.faces.vertices[iface + Side3D::B] = indexing::amr::sf(Side3D::B);
        cells.faces.vertices[iface + Side3D::T] = indexing::amr::sf(Side3D::T);
        cells.faces.vertices[iface + Side3D::Z] = indexing::amr::sf(Side3D::Z);
        cells.faces.vertices[iface + Side3D::F] = indexing::amr::sf(Side3D::F);

        for (index_t jn = 0; jn < n_nodes; ++jn) {
            cells.verts[ic * n_nodes + jn] = cube[jn];
        }
    }
    return cells;
}

} // namespace zephyr::geom::generator