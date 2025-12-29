#include <iostream>
#include <algorithm>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/utils/json.h>
#include <zephyr/mesh/side.h>

namespace zephyr::geom::generator {

using namespace mesh;

Cuboid::Cuboid(const Json& config)
    : Generator("cuboid"),
      m_xmin(0.0), m_xmax(1.0),
      m_ymin(0.0), m_ymax(1.0),
      m_zmin(0.0), m_zmax(1.0),
      m_nx(0), m_ny(0), m_nz(0) {

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
            ny = config["size"]["nz"].as<int>();
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
                std::string message = "Cuboid(json) error: Strange mesh sizes: " +
                                      std::to_string(nx) + " x " +
                                      std::to_string(ny) + " x " +
                                      std::to_string(nz) + "." +
                                      "Setup size.nx, size.ny, size.nz or each value";
                std::cerr << message << "\n";
                throw std::runtime_error(message);
            }
        }
    }
}

Cuboid::Cuboid(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) :
        Generator("cuboid"),
        m_xmin(xmin), m_xmax(xmax),
        m_ymin(ymin), m_ymax(ymax),
        m_zmin(zmin), m_zmax(zmax),
        m_nx(0), m_ny(0), m_nz(0), m_size(0),
        m_bounds() {
    check_params();
}

Box Cuboid::bbox() const {
    Vector3d vmin(m_xmin, m_ymin, m_zmin);
    Vector3d vmax(m_xmax, m_ymax, m_zmax);

    return {vmin, vmax};
}

void Cuboid::set_nx(int nx) {
    if (nx < 1) {
        throw std::runtime_error("Cuboid: Nx < 1");
    }
    m_nx = nx;
    m_ny = int(round(nx * (m_ymax - m_ymin) / (m_xmax - m_xmin)));
    m_nz = int(round(nx * (m_zmax - m_zmin) / (m_xmax - m_xmin)));
    compute_size();
}

void Cuboid::set_ny(int ny) {
    if (ny < 1) {
        throw std::runtime_error("Cuboid: Ny < 1");
    }
    m_ny = ny;
    m_nx = int(round(ny * (m_xmax - m_xmin) / (m_ymax - m_ymin)));
    m_nz = int(round(ny * (m_zmax - m_zmin) / (m_ymax - m_ymin)));
    compute_size();
}

void Cuboid::set_nz(int nz) {
    if (nz < 1) {
        throw std::runtime_error("Cuboid: Nz < 1");
    }
    m_nz = nz;
    m_nx = int(round(nz * (m_xmax - m_xmin) / (m_zmax - m_zmin)));
    m_ny = int(round(nz * (m_ymax - m_ymin) / (m_zmax - m_zmin)));
    compute_size();
}

void Cuboid::set_sizes(int nx, int ny, int nz) {
    if (nx < 1 || ny < 1 || nz < 1) {
        throw std::runtime_error("Cuboid: Nx < 1 or Ny < 1 or Nz < 1");
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
        std::string message = "Cuboid Error: x_min >= x_max";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (m_ymin >= m_ymax) {
        std::string message = "Cuboid Error: y_min >= y_max";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (m_zmin >= m_zmax) {
        std::string message = "Cuboid Error: zmin >= zmax";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    double dx = (m_xmax - m_xmin) / m_nx;
    double dy = (m_ymax - m_ymin) / m_ny;
    double dz = (m_ymax - m_ymin) / m_ny;

    double dmax = std::max(dx, std::max(dy, dz));
    double dmin = std::min(dx, std::min(dy, dz));
    if (dmax / dmin > 1.0e3) {
        std::cerr << "Cuboid Warning: Huge aspect ratio (> 1000)\n";
    }
}

void Cuboid::compute_size() {
    m_size = m_nx * m_ny * m_nz;
}

Grid Cuboid::make() {
    check_size(m_size);

    double dx = (m_xmax - m_xmin) / m_nx;
    double dy = (m_ymax - m_ymin) / m_ny;
    double dz = (m_zmax - m_zmin) / m_nz;

    Grid grid;
    grid.reserve_nodes((m_nx + 1) * (m_ny + 1) * (m_nz + 1));
    grid.reserve_cells(m_nx * m_ny * m_nz);

    std::vector nodes(m_nx + 1, std::vector(m_ny + 1, std::vector<NodeInput>(m_nz + 1)));
    for (int i = 0; i <= m_nx; ++i) {
        for (int j = 0; j <= m_ny; ++j) {
            for (int k = 0; k <= m_ny; ++k) {
                double x = m_xmin + i * dx;
                double y = m_ymin + j * dy;
                double z = m_zmin + k * dz;
                nodes[i][j][k] = Vector3d{x, y, z};

                grid.add_node(&nodes[i][j][k]);
            }
        }
    }

    std::vector<Boundary> bc(6);
    for (int i = 0; i < m_nx; ++i) {
        for (int j = 0; j < m_ny; ++j) {
            for (int k = 0; k < m_nz; ++k) {
                bc[Side3D::L] = i == 0 ? m_bounds.left :  Boundary::INNER;
                bc[Side3D::R] = i == m_nx-1 ? m_bounds.right : Boundary::INNER;
                bc[Side3D::B] = j == 0 ? m_bounds.bottom : Boundary::INNER;
                bc[Side3D::T] = j == m_ny-1 ? m_bounds.top : Boundary::INNER;
                bc[Side3D::Z] = k == 0 ? m_bounds.back : Boundary::INNER;
                bc[Side3D::F] = k == m_nz-1 ? m_bounds.front : Boundary::INNER;

                grid.add_cell(
                    CellType::HEXAHEDRON, {
                        &nodes[i][j][k],
                        &nodes[i + 1][j][k],
                        &nodes[i + 1][j + 1][k],
                        &nodes[i][j + 1][k],
                        &nodes[i][j][k+1],
                        &nodes[i + 1][j][k+1],
                        &nodes[i + 1][j + 1][k+1],
                        &nodes[i][j + 1][k+1]
                    }, bc);
            }
        }
    }
    return grid;
}

} // namespace zephyr::geom::generator