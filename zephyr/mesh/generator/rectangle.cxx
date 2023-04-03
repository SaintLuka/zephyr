#include <iostream>
#include <algorithm>

#include <zephyr/geom/cell.h>

#include <zephyr/geom/box.h>
#include <zephyr/mesh/generator/rectangle.h>


namespace zephyr { namespace mesh { namespace generator {

#ifdef ZEPHYR_ENABLE_MPI
using zephyr::network::mpi::Network;
#endif

#ifdef ZEPHYR_ENABLE_YAML
Rectangle::Rectangle(YAML::Node config)
    : Generator("rectangle"),
      m_xmin(0.0), m_xmax(1.0), m_ymin(0.0), m_ymax(1.0), m_nx(0), m_ny(0),
      m_left_flag(FaceFlag::WALL), m_right_flag(FaceFlag::WALL),
      m_bottom_flag(FaceFlag::WALL), m_top_flag(FaceFlag::WALL),
      m_voronoi(false) {

    if (!config["geometry"]) {
        throw std::runtime_error("Mesh config doesn't contain 'geometry'");
    }

    m_xmin = config["geometry"]["x_min"].as<double>();
    m_xmax = config["geometry"]["x_max"].as<double>();
    m_ymin = config["geometry"]["y_min"].as<double>();
    m_ymax = config["geometry"]["y_max"].as<double>();

    if (!config["boundary"]) {
        throw std::runtime_error("Mesh config doesn't contain 'boundary'");
    }
    m_left_flag = boundary_from_string(config["boundary"]["left"].as<std::string>());
    m_right_flag = boundary_from_string(config["boundary"]["right"].as<std::string>());
    m_bottom_flag = boundary_from_string(config["boundary"]["bottom"].as<std::string>());
    m_top_flag = boundary_from_string(config["boundary"]["top"].as<std::string>());

    if (config["voronoi"]) {
        m_voronoi = config["voronoi"].as<bool>();
    }

    if (config["cells"]) {
        set_size(config["cells"].as<int>());
    }
    else {
        int nx(0), ny(0);
        if (config["cells_per_x"]) {
            nx = config["cells_per_x"].as<int>();
        }
        if (config["cells_per_y"]) {
            ny = config["cells_per_y"].as<int>();
        }
        if (nx > 0 && ny > 0) {
            set_sizes(nx, ny);
        } else if (nx > 0) {
            set_nx(nx);
        } else if (ny > 0) {
            set_ny(ny);
        } else {
            throw std::runtime_error("Rectangle error: Strange mesh sizes");
        }
    }
}
#endif

Rectangle::Rectangle(double xmin, double xmax, double ymin, double ymax, bool voronoi) :
        Generator("rectangle"),
        m_xmin(xmin), m_xmax(xmax), m_ymin(ymin), m_ymax(ymax), m_nx(0), m_ny(0),
        m_left_flag(FaceFlag::UNDEFINED), m_right_flag(FaceFlag::UNDEFINED),
        m_bottom_flag(FaceFlag::UNDEFINED), m_top_flag(FaceFlag::UNDEFINED),
        m_voronoi(voronoi) {
    check_params();
}

Box Rectangle::bbox() const {
    Vector3d vmin(m_xmin, m_ymin, 0.0);
    Vector3d vmax(m_xmax, m_ymax, 0.0);

    return Box(vmin, vmax);
}

void Rectangle::set_nx(int nx) {
    if (nx < 1) {
        throw std::runtime_error("Rectangle: Nx < 1");
    }
    if (!m_voronoi) {
        m_nx = nx;
        m_ny = int(std::round(m_nx * (m_ymax - m_ymin) / (m_xmax - m_xmin)));
    }
    else {
        m_nx = int(std::round(nx * std::pow(3.0, 0.25)) / std::sqrt(2.0));
        m_ny = int(std::round(m_nx * (m_ymax - m_ymin) / (m_xmax - m_xmin) / std::sqrt(3.0)));
    }
    compute_size();
}

void Rectangle::set_ny(int ny) {
    if (ny < 1) {
        throw std::runtime_error("Rectangle: Ny < 1");
    }
    if (!m_voronoi) {
        m_ny = ny;
        m_nx = int(std::round(ny * (m_xmax - m_xmin) / (m_ymax - m_ymin)));
    }
    else {
        m_ny = int(std::round(ny / std::pow(3.0, 0.25)) / std::sqrt(2.0));
        m_nx = int(std::round(m_ny * (m_xmax - m_xmin) / (m_ymax - m_ymin) * std::sqrt(3.0)));
    }
    compute_size();
}

void Rectangle::set_sizes(int nx, int ny) {
    if (nx < 1 || ny < 1) {
        throw std::runtime_error("Rectangle: Nx < 1 or Ny < 1");
    }
    if (!m_voronoi) {
        m_nx = nx;
        m_ny = ny;
    }
    else {
        m_nx = int(std::round(nx * std::pow(3.0, 0.25)) / std::sqrt(2.0));
        m_ny = int(std::round(ny / std::pow(3.0, 0.25)) / std::sqrt(2.0));
    }
    compute_size();

    double dx = (m_xmax - m_xmin) / m_nx;
    double dy = (m_ymax - m_ymin) / m_ny;

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

    if (!m_voronoi) {
        double d = std::sqrt((m_xmax - m_xmin) * (m_ymax - m_ymin) / N);
        m_nx = int(std::round((m_xmax - m_xmin) / d));
        m_ny = int(std::round((m_ymax - m_ymin) / d));
    }
    else {
        double a = std::sqrt(2 * (m_xmax - m_xmin) * (m_ymax - m_ymin) / (std::sqrt(3.0) * N));
        double h = std::sqrt(3.0) * a;

        m_nx = int(std::round((m_xmax - m_xmin) / a));
        m_ny = int(std::round((m_ymax - m_ymin) / h));
    }
    compute_size();
    std::cout << m_nx << " " << m_ny << " " << m_size << "\n";
}

void Rectangle::set_boundary_flags(FaceFlag left, FaceFlag right, FaceFlag bottom, FaceFlag top) {
    m_left_flag = left;
    m_right_flag = right;
    m_bottom_flag = bottom;
    m_top_flag = top;
    if (periodic_along_x()) {
        m_left_flag = m_right_flag = FaceFlag::PERIODIC;
    }
    if (periodic_along_y()) {
        m_bottom_flag = m_top_flag = FaceFlag::PERIODIC;
    }
}

double Rectangle::x_min() const {
    return m_xmin;
}

double Rectangle::x_max() const {
    return m_xmax;
}

double Rectangle::y_min() const {
    return m_ymin;
}

double Rectangle::y_max() const {
    return m_ymax;
}

int Rectangle::nx() const {
    return m_nx;
}

int Rectangle::ny() const {
    return m_ny;
}

bool Rectangle::periodic_along_x() const {
    return m_left_flag == FaceFlag::PERIODIC || m_right_flag == FaceFlag::PERIODIC;
}

bool Rectangle::periodic_along_y() const {
    return m_bottom_flag == FaceFlag::PERIODIC || m_top_flag == FaceFlag::PERIODIC;
}

void Rectangle::check_params() const {
    if (m_xmin >= m_xmax) {
        std::cerr << "Rectangle Error: x_min >= x_max\n";
        throw std::runtime_error("Rectangle Error: x_min >= x_max");
    }
    if (m_ymin >= m_ymax) {
        std::cerr << "Rectangle Error: y_min >= y_max\n";
        throw std::runtime_error("Rectangle Error: y_min >= y_max");
    }
}

void Rectangle::compute_size() {
    if (!m_voronoi) {
        m_size = m_nx * m_ny;
    } else {
        m_size = 2 * m_nx * m_ny;
    }
    if (m_size > 1000000000) {
        std::cerr << "Attempt to create mesh with more than 1 billion cells\n";
        throw std::runtime_error("Attempt to create mesh with more than 1 billion cells");
    }
}

void Rectangle::initialize(Storage &cells) {
    cells.resize(m_size);

    if (m_voronoi) {
        init_voronoi(cells, Part(m_size));
    } else {
        init_classic(cells, Part(m_size));
    }
}

void Rectangle::initialize(Storage &cells, Part part) {
    if (m_voronoi) {
        init_voronoi(cells, part);
    } else {
        init_classic(cells, part);
    }
}

void Rectangle::init_classic(Storage &cells, Part part) const {
    using zephyr::geom::Side;
    using zephyr::geom::Cell;
    using zephyr::geom::ShortList2D;

    double dx = (m_xmax - m_xmin) / m_nx;
    double dy = (m_ymax - m_ymin) / m_ny;

    auto get_n = [this](int i, int j) -> int {
        return ((i + m_nx) % m_nx) * m_ny + (j + m_ny) % m_ny;
    };

    int counter = 0;
    for (int n = part.from; n < part.to; ++n, ++counter) {
        int i = n / m_ny;
        int j = n % m_ny;

        // Границы ячейки
        double x1 = m_xmin + i * dx;
        double x2 = m_xmin + (i + 1.0) * dx;

        double y1 = m_ymin + j * dy;
        double y2 = m_ymin + (j + 1.0) * dy;

        ShortList2D verts = {
                Vector3d(x1, y1, 0.0), Vector3d(x2, y1, 0.0),
                Vector3d(x1, y2, 0.0), Vector3d(x2, y2, 0.0)
        };

        Cell cell(verts);

        // Данные AMR
        cell.b_idx = n;
        cell.level   = 0;
        cell.flag    = 0;
        cell.z_idx       = 0;

        // Флаги граничных условий
        auto ordinary = FaceFlag::ORDINARY;

        cell.faces[Side::L].boundary = i > 0 ? ordinary : m_left_flag;
        cell.faces[Side::R].boundary = i < m_nx - 1 ? ordinary : m_right_flag;
        cell.faces[Side::B].boundary = j > 0 ? ordinary : m_bottom_flag;
        cell.faces[Side::T].boundary = j < m_ny - 1 ? ordinary : m_top_flag;

        cell.faces[Side::L].adjacent.rank = 0;
        cell.faces[Side::R].adjacent.rank = 0;
        cell.faces[Side::B].adjacent.rank = 0;
        cell.faces[Side::T].adjacent.rank = 0;

        cell.faces[Side::L].adjacent.ghost = 0;
        cell.faces[Side::R].adjacent.ghost = 0;
        cell.faces[Side::B].adjacent.ghost = 0;
        cell.faces[Side::T].adjacent.ghost = 0;

        cell.faces[Side::L].adjacent.index = get_n(i - 1, j);
        cell.faces[Side::R].adjacent.index = get_n(i + 1, j);
        cell.faces[Side::B].adjacent.index = get_n(i, j - 1);
        cell.faces[Side::T].adjacent.index = get_n(i, j + 1);

        cells[counter].geom() = cell;
    }
}

// Центр описаной окружности
Vector3d circle_center(const Vector3d& A, const Vector3d& B, const Vector3d& C) {
    double D = 2 * (A.x() * (B.y() - C.y()) + B.x() * (C.y() - A.y()) + C.x() * (A.y() - B.y()));
    double Ux = (A.squaredNorm() * (B.y() - C.y()) +
                 B.squaredNorm() * (C.y() - A.y()) +
                 C.squaredNorm() * (A.y() - B.y())) / D;
    double Uy = (A.squaredNorm() * (C.x() - B.x()) +
                 B.squaredNorm() * (A.x() - C.x()) +
                 C.squaredNorm() * (B.x() - A.x())) / D;
    return Vector3d(Ux, Uy, 0.0);
}

void Rectangle::init_voronoi(Storage &cells, Part part) const {
    using zephyr::geom::Side;
    using zephyr::geom::Cell;
    using zephyr::geom::ShortList2D;

    double Lx = m_xmax - m_xmin;
    double Ly = m_ymax - m_ymin;

    double a = Lx / m_nx;
    double h = Ly / m_ny;

    auto get_n = [this](int i, int j, int s) -> int {
        return 2 * (((i + m_nx) % m_nx) * m_ny + (j + m_ny) % m_ny) + s;
    };

    auto get_i = [this](int n) -> int { return (n / 2) / m_ny; };
    auto get_j = [this](int n) -> int { return (n / 2) % m_ny; };
    auto get_s = [this](int n) -> int { return n % 2; };
    auto inv = [](int s) -> int { return (s + 1) % 2; };

    // Центры ячеек
    auto center = [this, a, h](int i, int j, int s) -> Vector3d {
        i = (i + m_nx) % m_nx;
        j = (j + m_ny) % m_ny;
        return {
            m_xmin + (i + 0.5 * s + 0.25) * a,
            m_ymin + (j + 0.5 * s + 0.25) * h,
            0.0 };
    };

    int counter = 0;
    for (int n = part.from; n < part.to; ++n, ++counter) {
        int i = get_i(n);
        int j = get_j(n);
        int s = get_s(n);

        // Центр ячейки
        Vector3d c = center(i, j, s);

        std::vector<int> neibs;

        std::vector<int> neibs_n;
        neibs_n.reserve(6);
        std::vector<Vector3d> neibs_c;
        neibs_c.reserve(6);

        Vector3d nc;

        // Сосед слева
        neibs_n.push_back(get_n(i - 1, j, s));
        nc = center(i - 1, j, s);
        if (nc.x() > c.x()) {
            nc.x() -= Lx;
        }
        neibs_c.push_back(nc);

        // Сосед внизу слева
        neibs_n.push_back(get_n(i - inv(s), j - inv(s), inv(s)));
        nc = center(i - inv(s), j - inv(s), inv(s));
        if (nc.x() > c.x()) { nc.x() -= Lx; }
        if (nc.y() > c.y()) { nc.y() -= Ly; }
        neibs_c.push_back(nc);

        // Сосед внизу справа
        neibs_n.push_back(get_n(i + s, j - inv(s), inv(s)));
        nc = center(i + s, j - inv(s), inv(s));
        if (nc.x() < c.x()) { nc.x() += Lx; }
        if (nc.y() > c.y()) { nc.y() -= Ly; }
        neibs_c.push_back(nc);

        // Сосед справа
        neibs_n.push_back(get_n(i + 1, j, s));
        nc = center(i + 1, j, s);
        if (nc.x() < c.x()) { nc.x() += Lx; }
        neibs_c.push_back(nc);

        // Сосед сверху справа
        neibs_n.push_back(get_n(i + s, j + s, inv(s)));
        nc = center(i + s, j + s, inv(s));
        if (nc.x() < c.x()) { nc.x() += Lx; }
        if (nc.y() < c.y()) { nc.y() += Ly; }
        neibs_c.push_back(nc);

        // Сосед сверху слева
        neibs_n.push_back(get_n(i - inv(s), j + s, inv(s)));
        nc = center(i - inv(s), j + s, inv(s));
        if (nc.x() > c.x()) { nc.x() -= Lx; }
        if (nc.y() < c.y()) { nc.y() += Ly; }
        neibs_c.push_back(nc);

        std::vector<Vector3d> verts(6);
        for (int k1 = 0; k1 < 6; ++k1) {
            int k2 = (k1 + 5) % 6;

            verts[k1] = circle_center(c, neibs_c[k1], neibs_c[k2]);
        }

        Cell cell(verts);

        // Данные AMR
        cell.b_idx = n;
        cell.level   = 0;
        cell.flag    = 0;
        cell.z_idx       = 0;

        // Флаги граничных условий
        auto ordinary = FaceFlag::ORDINARY;

        cell.faces[Side::L].boundary = ordinary;
        cell.faces[Side::R].boundary = ordinary;
        cell.faces[Side::B].boundary = ordinary;
        cell.faces[Side::T].boundary = ordinary;
        cell.faces[Side::X].boundary = ordinary;
        cell.faces[Side::F].boundary = ordinary;

        cell.faces[Side::L].adjacent.rank = 0;
        cell.faces[Side::R].adjacent.rank = 0;
        cell.faces[Side::B].adjacent.rank = 0;
        cell.faces[Side::T].adjacent.rank = 0;
        cell.faces[Side::X].adjacent.rank = 0;
        cell.faces[Side::F].adjacent.rank = 0;

        cell.faces[Side::L].adjacent.ghost = 0;
        cell.faces[Side::R].adjacent.ghost = 0;
        cell.faces[Side::B].adjacent.ghost = 0;
        cell.faces[Side::T].adjacent.ghost = 0;
        cell.faces[Side::X].adjacent.ghost = 0;
        cell.faces[Side::F].adjacent.ghost = 0;

        cell.faces[Side::L].adjacent.index = neibs_n[0];
        cell.faces[Side::R].adjacent.index = neibs_n[1];
        cell.faces[Side::B].adjacent.index = neibs_n[2];
        cell.faces[Side::T].adjacent.index = neibs_n[3];
        cell.faces[Side::X].adjacent.index = neibs_n[4];
        cell.faces[Side::F].adjacent.index = neibs_n[5];

        cells[counter].geom() = cell;
    }
}

} // generator
} // mesh
} // zephyr