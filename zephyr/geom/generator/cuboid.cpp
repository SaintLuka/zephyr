#include <iostream>
#include <algorithm>

#include <zephyr/geom/box.h>
#include <zephyr/geom/primitives/amr_cell.h>

#include <zephyr/geom/generator/cuboid.h>

namespace zephyr::geom::generator {

#ifdef ZEPHYR_ENABLE_YAML
Cuboid::Cuboid(YAML::Node config)
    : Generator("cuboid"),
      m_xmin(0.0), m_xmax(1.0),
      m_ymin(0.0), m_ymax(1.0),
      m_zmin(0.0), m_zmax(1.0),
      m_nx(0), m_ny(0), m_nz(0),
      m_left_flag(Boundary::UNDEFINED), m_right_flag(Boundary::UNDEFINED),
      m_bottom_flag(Boundary::UNDEFINED), m_top_flag(Boundary::UNDEFINED),
      m_back_flag(Boundary::UNDEFINED), m_front_flag(Boundary::UNDEFINED) {

    if (!config["geometry"]) {
        throw std::runtime_error("Mesh config doesn't contain 'geometry'");
    }

    m_xmin = config["geometry"]["x_min"].as<double>();
    m_xmax = config["geometry"]["x_max"].as<double>();
    m_ymin = config["geometry"]["y_min"].as<double>();
    m_ymax = config["geometry"]["y_max"].as<double>();
    m_zmin = config["geometry"]["z_min"].as<double>();
    m_zmax = config["geometry"]["z_max"].as<double>();

    if (!config["boundary"]) {
        throw std::runtime_error("Mesh config doesn't contain 'boundary'");
    }
    m_left_flag = boundary_from_string(config["boundary"]["left"].as<std::string>());
    m_right_flag = boundary_from_string(config["boundary"]["right"].as<std::string>());
    m_bottom_flag = boundary_from_string(config["boundary"]["bottom"].as<std::string>());
    m_top_flag = boundary_from_string(config["boundary"]["top"].as<std::string>());
    m_back_flag = boundary_from_string(config["boundary"]["back"].as<std::string>());
    m_front_flag = boundary_from_string(config["boundary"]["front"].as<std::string>());

    if (config["cells"]) {
        set_size(config["cells"].as<size_t>());
    }
    else {
        size_t nx(0), ny(0), nz(0);
        if (config["cells_per_x"]) {
            nx = config["cells_per_x"].as<size_t>();
        }
        if (config["cells_per_y"]) {
            ny = config["cells_per_y"].as<size_t>();
        }
        if (config["cells_per_z"]) {
            ny = config["cells_per_z"].as<size_t>();
        }
        if (nx > 0 && ny > 0 && nz > 0) {
            set_sizes(nx, ny, nz);
        }
        else {
            if (nx > 0) {
                set_nx(nx);
            }
            else if (ny > 0) {
                set_ny(ny);
            }
            else if (nz > 0) {
                set_nz(nz);
            }
            else {
                throw std::runtime_error("Cuboid error: Strange mesh sizes");
            }
        }
    }
}
#endif

Cuboid::Cuboid(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) :
        Generator("cuboid"),
        m_xmin(xmin), m_xmax(xmax),
        m_ymin(ymin), m_ymax(ymax),
        m_zmin(zmin), m_zmax(zmax),
        m_nx(0), m_ny(0), m_nz(0),
        m_left_flag(Boundary::UNDEFINED), m_right_flag(Boundary::UNDEFINED),
        m_bottom_flag(Boundary::UNDEFINED), m_top_flag(Boundary::UNDEFINED),
        m_back_flag(Boundary::UNDEFINED), m_front_flag(Boundary::UNDEFINED) {
    check_params();
}

Box Cuboid::bbox() const {
    Vector3d vmin(m_xmin, m_ymin, m_zmin);
    Vector3d vmax(m_xmax, m_ymax, m_zmax);

    return Box(vmin, vmax);
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

void Cuboid::set_boundary_flags(
        Boundary left, Boundary right,
        Boundary bottom, Boundary top,
        Boundary back, Boundary front) {

    m_left_flag = left;
    m_right_flag = right;
    m_bottom_flag = bottom;
    m_top_flag = top;
    m_back_flag = back;
    m_front_flag = front;
    if (periodic_along_x()) {
        m_left_flag = m_right_flag = Boundary::PERIODIC;
    }
    if (periodic_along_y()) {
        m_bottom_flag = m_top_flag = Boundary::PERIODIC;
    }
    if (periodic_along_z()) {
        m_back_flag = m_front_flag = Boundary::PERIODIC;
    }
}

double Cuboid::xmin() const {
    return m_xmin;
}

double Cuboid::xmax() const {
    return m_xmax;
}

double Cuboid::ymin() const {
    return m_ymin;
}

double Cuboid::ymax() const {
    return m_ymax;
}

double Cuboid::zmin() const {
    return m_ymin;
}

double Cuboid::zmax() const {
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

bool Cuboid::periodic_along_x() const {
    return m_left_flag == Boundary::PERIODIC || m_right_flag == Boundary::PERIODIC;
}

bool Cuboid::periodic_along_y() const {
    return m_bottom_flag == Boundary::PERIODIC || m_top_flag == Boundary::PERIODIC;
}

bool Cuboid::periodic_along_z() const {
    return m_back_flag == Boundary::PERIODIC || m_front_flag == Boundary::PERIODIC;
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

Grid Cuboid::create() {
    /*
    using zephyr::geom::Side;
    using zephyr::geom::AmrCell;
    using zephyr::geom::ShortList3D;

    double dx = (m_xmax - m_xmin) / m_nx;
    double dy = (m_ymax - m_ymin) / m_ny;
    double dz = (m_zmax - m_zmin) / m_nz;

    cells.resize(m_size);

    for (int n = 0; n < cells.size(); ++n) {
        int i = (n / m_nz) / m_ny;
        int j = (n / m_nz) % m_ny;
        int k = n % m_nz;

        // Границы ячейки
        double x1 = m_xmin + i * dx;
        double x2 = m_xmin + (i + 1.0) * dx;

        double y1 = m_ymin + j * dy;
        double y2 = m_ymin + (j + 1.0) * dy;

        double z1 = m_zmin + k * dz;
        double z2 = m_zmin + (k + 1.0) * dz;

        ShortList3D verts = {
                Vector3d(x1, y1, z1), Vector3d(x2, y1, z1),
                Vector3d(x1, y2, z1), Vector3d(x2, y2, z1),
                Vector3d(x1, y1, z2), Vector3d(x2, y1, z2),
                Vector3d(x1, y2, z2), Vector3d(x2, y2, z2)
        };

        AmrCell g_cell(verts);

        auto ordinary = Boundary::ORDINARY;

        g_cell.faces[Side::L].boundary = i > 0 ? ordinary : m_left_flag;
        g_cell.faces[Side::R].boundary = i < m_nx - 1 ? ordinary : m_right_flag;
        g_cell.faces[Side::B].boundary = j > 0 ? ordinary : m_bottom_flag;
        g_cell.faces[Side::T].boundary = j < m_ny - 1 ? ordinary : m_top_flag;
        g_cell.faces[Side::X].boundary = k > 0 ? ordinary : m_back_flag;
        g_cell.faces[Side::F].boundary = k < m_nz - 1 ? ordinary : m_front_flag;

        cells[n].geom() = g_cell;
    }
     */
}

} // namespace zephyr::geom::generator