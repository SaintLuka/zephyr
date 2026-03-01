#include <iostream>
#include <format>
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
      m_xmin(0.0), m_xmax(1.0),
      m_ymin(0.0), m_ymax(1.0),
      m_zmin(0.0), m_zmax(1.0) {

    if (!config["geometry"]) {
        throw std::runtime_error("Cuboid config doesn't contain key 'geometry'");
    }

    m_xmin = config["geometry"]["x_min"].as<double>();
    m_xmax = config["geometry"]["x_max"].as<double>();
    m_ymin = config["geometry"]["y_min"].as<double>();
    m_ymax = config["geometry"]["y_max"].as<double>();
    m_zmin = config["geometry"]["z_min"].as<double>();
    m_zmax = config["geometry"]["z_max"].as<double>();

    if (!config["bounds"]) {
        throw std::runtime_error("Cuboid config doesn't contain key 'bounds'");
    }
    m_bounds.left   = boundary_from_string(config["bounds"]["left"].as<std::string>());
    m_bounds.right  = boundary_from_string(config["bounds"]["right"].as<std::string>());
    m_bounds.bottom = boundary_from_string(config["bounds"]["bottom"].as<std::string>());
    m_bounds.top    = boundary_from_string(config["bounds"]["top"].as<std::string>());
    m_bounds.back   = boundary_from_string(config["bounds"]["back"].as<std::string>());
    m_bounds.front  = boundary_from_string(config["bounds"]["front"].as<std::string>());

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

Cuboid::Cuboid(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) :
        Generator("cuboid"),
        m_xmin(xmin), m_xmax(xmax),
        m_ymin(ymin), m_ymax(ymax),
        m_zmin(zmin), m_zmax(zmax) {
    check_params();
}

Box Cuboid::bbox() const {
    Vector3d vmin(m_xmin, m_ymin, m_zmin);
    Vector3d vmax(m_xmax, m_ymax, m_zmax);

    return {vmin, vmax};
}

void Cuboid::set_nx(int nx) {
    if (nx < 1) {
        throw std::runtime_error("Cuboid::set_nx: Nx < 1");
    }
    m_nx = nx;
    m_ny = int(round(nx * (m_ymax - m_ymin) / (m_xmax - m_xmin)));
    m_nz = int(round(nx * (m_zmax - m_zmin) / (m_xmax - m_xmin)));
    compute_size();
}

void Cuboid::set_ny(int ny) {
    if (ny < 1) {
        throw std::runtime_error("Cuboid::set_ny: Ny < 1");
    }
    m_ny = ny;
    m_nx = int(round(ny * (m_xmax - m_xmin) / (m_ymax - m_ymin)));
    m_nz = int(round(ny * (m_zmax - m_zmin) / (m_ymax - m_ymin)));
    compute_size();
}

void Cuboid::set_nz(int nz) {
    if (nz < 1) {
        throw std::runtime_error("Cuboid::set_nz: Nz < 1");
    }
    m_nz = nz;
    m_nx = int(round(nz * (m_xmax - m_xmin) / (m_zmax - m_zmin)));
    m_ny = int(round(nz * (m_ymax - m_ymin) / (m_zmax - m_zmin)));
    compute_size();
}

void Cuboid::set_sizes(int nx, int ny, int nz) {
    if (nx < 1 || ny < 1 || nz < 1) {
        throw std::runtime_error("Cuboid::set_sizes: Nx < 1 or Ny < 1 or Nz < 1");
    }
    m_nx = nx;
    m_ny = ny;
    m_nz = nz;
    compute_size();

    double dx = (m_xmax - m_xmin) / m_nx;
    double dy = (m_ymax - m_ymin) / m_ny;
    double dz = (m_zmax - m_zmin) / m_nz;

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
    double d = std::cbrt((m_xmax - m_xmin) * (m_ymax - m_ymin) * (m_zmax - m_zmin) / N);
    m_nx = int(round((m_xmax - m_xmin) / d));
    m_ny = int(round((m_ymax - m_ymin) / d));
    m_nz = int(round((m_zmax - m_zmin) / d));
    compute_size();
}

void Cuboid::set_boundaries(Boundaries bounds) {
    m_bounds = bounds;

    if (periodic_along_x()) {
        m_bounds.left = m_bounds.right = Boundary::PERIODIC;
    }
    if (periodic_along_y()) {
        m_bounds.bottom = m_bounds.top = Boundary::PERIODIC;
    }
    if (periodic_along_z()) {
        m_bounds.back = m_bounds.front = Boundary::PERIODIC;
    }
}

double Cuboid::x_min() const {
    return m_xmin;
}

double Cuboid::x_max() const {
    return m_xmax;
}

double Cuboid::y_min() const {
    return m_ymin;
}

double Cuboid::y_max() const {
    return m_ymax;
}

double Cuboid::z_min() const {
    return m_ymin;
}

double Cuboid::z_max() const {
    return m_ymax;
}

int Cuboid::nx() const {
    return m_nx;
}

int Cuboid::ny() const {
    return m_ny;
}

int Cuboid::nz() const {
    return m_ny;
}

Cuboid::Boundaries Cuboid::bounds() const {
    return m_bounds;
}

bool Cuboid::periodic_along_x() const {
    return m_bounds.left == Boundary::PERIODIC || m_bounds.right == Boundary::PERIODIC;
}

bool Cuboid::periodic_along_y() const {
    return m_bounds.bottom == Boundary::PERIODIC || m_bounds.top == Boundary::PERIODIC;
}

bool Cuboid::periodic_along_z() const {
    return m_bounds.back == Boundary::PERIODIC || m_bounds.front == Boundary::PERIODIC;
}

void Cuboid::check_params() const {
    if (m_xmin >= m_xmax) {
        throw std::runtime_error("Cuboid::check_params: x_min >= x_max");
    }
    if (m_ymin >= m_ymax) {
        throw std::runtime_error("Cuboid::check_params: y_min >= y_max");
    }
    if (m_zmin >= m_zmax) {
        throw std::runtime_error("Cuboid::check_params: z_min >= z_max");
    }
    double dx = (m_xmax - m_xmin) / m_nx;
    double dy = (m_ymax - m_ymin) / m_ny;
    double dz = (m_ymax - m_ymin) / m_ny;

    double dmax = std::max(dx, std::max(dy, dz));
    double dmin = std::min(dx, std::min(dy, dz));
    if (dmax / dmin > 1.0e3) {
        std::cerr << "Cuboid::check_params warning: Huge aspect ratio (> 1000)\n";
    }
}

void Cuboid::compute_size() {
    m_size = m_nx * m_ny * m_nz;
}

Grid Cuboid::make() const {
    check_size(m_size);

    double dx = (m_xmax - m_xmin) / m_nx;
    double dy = (m_ymax - m_ymin) / m_ny;
    double dz = (m_zmax - m_zmin) / m_nz;

    Grid grid;
    grid.reserve_nodes((m_nx + 1) * (m_ny + 1) * (m_nz + 1));
    grid.reserve_cells(m_nx * m_ny * m_nz);

    std::vector nodes(m_nx + 1, std::vector(m_ny + 1, std::vector<GridNode::Ptr>(m_nz + 1)));
    for (int i = 0; i <= m_nx; ++i) {
        for (int j = 0; j <= m_ny; ++j) {
            for (int k = 0; k <= m_nz; ++k) {
                double x = m_xmin + i * dx;
                double y = m_ymin + j * dy;
                double z = m_zmin + k * dz;
                nodes[i][j][k] = GridNode::create({x, y, z});
            }
        }
    }

    std::vector<Boundary> bc(6);
    std::vector<GridNode::Ptr> cube_nodes(8);
    for (int i = 0; i < m_nx; ++i) {
        bc[Side3D::L] = i == 0 ?      m_bounds.left :  Boundary::INNER;
        bc[Side3D::R] = i == m_nx-1 ? m_bounds.right : Boundary::INNER;
        for (int j = 0; j < m_ny; ++j) {
            bc[Side3D::B] = j == 0 ?      m_bounds.bottom : Boundary::INNER;
            bc[Side3D::T] = j == m_ny-1 ? m_bounds.top : Boundary::INNER;
            for (int k = 0; k < m_nz; ++k) {
                bc[Side3D::Z] = k == 0 ? m_bounds.back : Boundary::INNER;
                bc[Side3D::F] = k == m_nz-1 ? m_bounds.front : Boundary::INNER;

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

void Cuboid::initialize(AmrCells& cells) const {
    bool x_period = periodic_along_x();
    bool y_period = periodic_along_y();
    bool z_period = periodic_along_z();

    double hx = (m_xmax - m_xmin) / m_nx;
    double hy = (m_ymax - m_ymin) / m_ny;
    double hz = (m_zmax - m_zmin) / m_nz;

    auto get_index = [=, this](index_t i, index_t j, index_t k) -> index_t {
        return m_nz * (m_ny * i + j) + k;
    };

    auto get_index_pair = [=, this](index_t n) -> std::array<index_t, 3> {
        return {(n / m_nz) / m_ny, (n / m_nz) % m_ny, n % m_nz};
    };

    auto get_vertex = [=, this](index_t i, index_t j, index_t k) -> Vector3d {
        return {
                m_xmin + ((m_xmax - m_xmin) * i) / m_nx,
                m_ymin + ((m_ymax - m_ymin) * j) / m_ny,
                m_zmin + ((m_zmax - m_zmin) * k) / m_nz
        };
    };

    auto neib_index = [=, this](index_t i, index_t j, index_t k, Side3D side) -> index_t {
        if (side == Side3D::LEFT) {
            return i == 0 && !x_period ?  get_index(i, j, k) : get_index((i - 1 + m_nx) % m_nx, j, k);
        }
        if (side == Side3D::RIGHT) {
            return i == m_nx - 1 && !x_period ? get_index(i, j, k) : get_index((i + 1) % m_nx, j, k);
        }
        if (side == Side3D::BOTTOM) {
            return j == 0 && !y_period ? get_index(i, j, k) : get_index(i, (j - 1 + m_ny) % m_ny, k);
        }
        if (side == Side3D::TOP) {
            return j == m_ny - 1 && !y_period ? get_index(i, j, k): get_index(i, (j + 1) % m_ny, k);
        }
        if (side == Side3D::BACK) {
            return k == 0 && !z_period ? get_index(i, j, k) : get_index(i, j, (k - 1 + m_nz) % m_nz);
        }
        if (side == Side3D::FRONT) {
            return k == m_nz - 1 && !z_period ? get_index(i, j, k): get_index(i, j, (k + 1) % m_nz);
        }
        throw std::runtime_error("Strange side #265");
    };

    cells.set_dimension(3);
    cells.set_adaptive(true);
    cells.set_linear(true);
    cells.set_axial(false);

    cells.resize_amr(m_size);

    int n_faces = 24;
    int n_nodes = 27;

    for (index_t ic = 0; ic < m_size; ++ic) {
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
        cells.face_begin[ic] = ic * n_faces;
        cells.node_begin[ic] = ic * n_nodes;
        cells.face_begin[ic + 1] = (ic + 1) * n_faces;
        cells.node_begin[ic + 1] = (ic + 1) * n_nodes;

        // INIT FACES
        for (auto iface: cells.faces_range(ic)) {
            cells.faces.set_undefined(iface);
            cells.faces.area_alt[iface] = NAN;
        }

        index_t iface = ic * n_faces;

        cells.faces.boundary[iface + Side3D::L] = i > 0 ? Boundary::INNER : m_bounds.left;
        cells.faces.boundary[iface + Side3D::R] = i < m_nx - 1 ? Boundary::INNER : m_bounds.right;
        cells.faces.boundary[iface + Side3D::B] = j > 0 ? Boundary::INNER : m_bounds.bottom;
        cells.faces.boundary[iface + Side3D::T] = j < m_ny - 1 ? Boundary::INNER : m_bounds.top;
        cells.faces.boundary[iface + Side3D::Z] = k > 0 ? Boundary::INNER : m_bounds.back;
        cells.faces.boundary[iface + Side3D::F] = k < m_nz - 1 ? Boundary::INNER : m_bounds.front;

        for (auto side: Side3D::items()) {
            cells.faces.adjacent.rank[iface + side] = 0;
            cells.faces.adjacent.index[iface + side] = neib_index(i, j, k, side);
            cells.faces.adjacent.alien[iface + side] = -1;
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
}

} // namespace zephyr::geom::generator