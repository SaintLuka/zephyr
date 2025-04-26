#include <iostream>
#include <algorithm>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/utils/json.h>

namespace zephyr::geom::generator {

Rectangle::Rectangle(const Json& config)
    : Generator("rectangle"),
      m_xmin(0.0), m_xmax(1.0), m_ymin(0.0), m_ymax(1.0), m_nx(0), m_ny(0),
      m_voronoi(false) {

    m_axial = false;
    if (config["axial"]) {
        m_axial = config["axial"].as<bool>();
    }

    if (!config["geometry"]) {
        throw std::runtime_error("Rectangle config doesn't contain key 'geometry'");
    }

    m_xmin = config["geometry"]["x_min"].as<double>();
    m_xmax = config["geometry"]["x_max"].as<double>();
    m_ymin = config["geometry"]["y_min"].as<double>();
    m_ymax = config["geometry"]["y_max"].as<double>();

    if (!config["bounds"]) {
        throw std::runtime_error("Rectangle config doesn't contain key 'bounds'");
    }
    m_bounds.left   = boundary_from_string(config["bounds"]["left"].as<std::string>());
    m_bounds.right  = boundary_from_string(config["bounds"]["right"].as<std::string>());
    m_bounds.bottom = boundary_from_string(config["bounds"]["bottom"].as<std::string>());
    m_bounds.top    = boundary_from_string(config["bounds"]["top"].as<std::string>());

    if (config["voronoi"]) {
        m_voronoi = config["voronoi"].as<bool>();
    }

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
            std::string message = "Rectangle(json) error: Strange mesh sizes: " +
                                  std::to_string(nx) + " x " + std::to_string(ny) + "." +
                                  "Setup size.nx, size.ny or both";
            std::cerr << message << "\n";
            throw std::runtime_error(message);
        }
    }
}

Rectangle::Rectangle()
    : Rectangle(0.0, 1.0, 0.0, 1.0, false) {
    set_boundaries(Boundaries{});
}

Rectangle::Rectangle(double xmin, double xmax, double ymin, double ymax, bool voronoi) :
        Generator("rectangle"),
        m_nx(0),
        m_ny(0),
        m_size(0),
        m_xmin(xmin),
        m_xmax(xmax),
        m_ymin(ymin),
        m_ymax(ymax),
        m_bounds(),
        m_voronoi(voronoi) {
    check_params();
}

int Rectangle::size() const {
    return m_size;
}

Box Rectangle::bbox() const {
    Vector3d vmin(m_xmin, m_ymin, 0.0);
    Vector3d vmax(m_xmax, m_ymax, 0.0);

    return {vmin, vmax};
}

void Rectangle::set_axial(bool axial) {
    m_axial = axial;
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
}

void Rectangle::set_boundaries(Boundaries bounds) {
    m_bounds = bounds;
    if (periodic_along_x()) {
        m_bounds.left = m_bounds.right = Boundary::PERIODIC;
    }
    if (periodic_along_y()) {
        m_bounds.bottom = m_bounds.top = Boundary::PERIODIC;
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

Rectangle::Boundaries Rectangle::bounds() const {
    return m_bounds;
}

bool Rectangle::periodic_along_x() const {
    return m_bounds.left == Boundary::PERIODIC || m_bounds.right == Boundary::PERIODIC;
}

bool Rectangle::periodic_along_y() const {
    return m_bounds.bottom == Boundary::PERIODIC || m_bounds.top == Boundary::PERIODIC;
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

Grid Rectangle::make() {
    check_size();

    if (m_axial && m_ymin == 0.0) {
        m_bounds.bottom = Boundary::WALL;
    }

    if (m_voronoi) {
        return create_voronoi();
    } else {
        return create_classic();
    }
}

Grid Rectangle::create_classic() const {
    double dx = (m_xmax - m_xmin) / m_nx;
    double dy = (m_ymax - m_ymin) / m_ny;

    std::vector<std::vector<GNode::Ptr>> nodes(
            m_nx + 1,
            std::vector<GNode::Ptr>(m_ny + 1, nullptr));

    int n_nodes = 0;
    for (int i = 0; i <= m_nx; ++i) {
        for (int j = 0; j <= m_ny; ++j) {
            double x = m_xmin + i * dx;
            double y = m_ymin + j * dy;
            nodes[i][j] = GNode::create(x, y);
            nodes[i][j]->index = n_nodes;
            ++n_nodes;
        }
    }

    for (int i = 0; i <= m_nx; ++i) {
        nodes[i][0]->add_boundary(m_bounds.bottom);
        nodes[i][m_ny]->add_boundary(m_bounds.top);
    }
    for (int j = 0; j <= m_ny; ++j) {
        nodes[0][j]->add_boundary(m_bounds.left);
        nodes[m_nx][j]->add_boundary(m_bounds.right);
    }

    Grid grid;
    grid.set_axial(m_axial);

    grid.reserve_nodes((m_nx + 1) * (m_ny + 1));
    for (int i = 0; i <= m_nx; ++i) {
        for (int j = 0; j <= m_ny; ++j) {
            grid += nodes[i][j];
        }
    }

    grid.reserve_cells(m_nx * m_ny);
    int n_cells = 0;
    for (int i = 0; i < m_nx; ++i) {
        for (int j = 0; j < m_ny; ++j) {
            GCell cell = GCell::quad(
                    {
                            nodes[i][j],
                            nodes[i + 1][j],
                            nodes[i + 1][j + 1],
                            nodes[i][j + 1]
                    });
            cell.index = n_cells;
            ++n_cells;

            grid += cell;
        }
    }

    grid.setup_adjacency();
    grid.assume_structured(m_nx, m_ny);

    return grid;
}

Grid Rectangle::create_voronoi() const {
    size_t Nb = m_nx * m_ny;
    double DX = m_xmax - m_xmin;
    double DY = m_ymax - m_ymin;
    size_t Ny = size_t(std::floor(std::sqrt(std::sqrt(3.0) * DY * Nb / 2.0 / DX))) + 1;
    double h = DY / 2.0 / Ny;
    double D = h / std::sqrt(0.75);

    size_t Nx = size_t(std::floor(DX / 1.5 / D));

    double x_shift = m_xmin + (DX - Nx * 1.5 * D) / 2.0 - 0.5 * D;
    double y_shift = m_ymin;

    // Вершины в виде таблицы
    std::vector<std::vector<GNode::Ptr>> vertices(
            Nx + 2, std::vector<GNode::Ptr>(2 * Ny + 1, nullptr)
    );

    for (size_t j = 0; j <= 2 * Ny; ++j) {
        double y = y_shift + h * j;

        // Часть вершин на левой границе пропускаем
        if (j % 2 == 0) {
            vertices[0][j] = GNode::create(m_xmin, y);
            vertices[0][j]->add_boundary(m_bounds.left);
        }

        for (size_t i = 1; i <= Nx; ++i) {
            double x = x_shift + double(3 * i - (i + j) % 2) * 0.5 * D;
            vertices[i][j] = GNode::create(x, y);

            if (j == 0) {
                vertices[i][j]->add_boundary(m_bounds.bottom);
            } else if (j == 2 * Ny) {
                vertices[i][j]->add_boundary(m_bounds.top);
            }
        }

        // Часть вершин на правой границе пропускаем
        if (j == 0 || j == 2 * Ny || j % 2 == Nx % 2) {
            vertices[Nx + 1][j] = GNode::create(m_xmax, y);
            vertices[Nx + 1][j]->add_boundary(m_bounds.right);
        }
    }

    Grid grid;

    // Переносим вершины в один массив
    grid.reserve_nodes((Nx + 2)*(2*Ny+1));
    int n_nodes = 0;
    for (size_t i = 0; i <= Nx + 1; ++i) {
        for (size_t j = 0; j <= 2 * Ny; ++j) {
            if (vertices[i][j]) {
                vertices[i][j]->index = n_nodes;
                ++n_nodes;

                grid += vertices[i][j];
            }
        }
    }

    using VList = std::vector<GNode::Ptr>;

    auto erase_nulls = [](VList& vlist) {
        auto to_remove = std::remove(vlist.begin(), vlist.end(), nullptr);
        vlist.erase(to_remove, vlist.end());
    };

    // Создаем массив ячеек
    grid.reserve_nodes(m_size);
    int n_cells = 0;
    for (size_t i = 0; i <= Nx; ++i) {
        if (i % 2 == 0) {
            for (size_t j = 0; j < 2 * Ny - 1; j += 2) {
                VList vlist = {
                        vertices[i][j], vertices[i + 1][j],
                        vertices[i + 1][j + 1], vertices[i + 1][j + 2],
                        vertices[i][j + 2], vertices[i][j + 1]
                };
                erase_nulls(vlist);
                GCell cell = GCell::polygon(vlist);
                cell.index = n_cells;
                grid += cell;
                ++n_cells;
            }
        } else {
            // j == 0
            auto cell = GCell::quad(
                    {
                            vertices[i][0], vertices[i + 1][0],
                            vertices[i + 1][1], vertices[i][1]
                    });
            cell.index = n_cells;
            grid += cell;
            ++n_cells;

            for (size_t j = 1; j < 2 * Ny - 2; j += 2) {
                VList vlist = {
                        vertices[i][j], vertices[i + 1][j],
                        vertices[i + 1][j + 1], vertices[i + 1][j + 2],
                        vertices[i][j + 2], vertices[i][j + 1]
                };
                erase_nulls(vlist);
                cell = GCell::polygon(vlist);
                cell.index = n_cells;
                grid += cell;
                ++n_cells;
            }

            // j == last; i % 2 == 1
            cell = GCell::quad(
                    {
                            vertices[i][2 * Ny - 1], vertices[i + 1][2 * Ny - 1],
                            vertices[i + 1][2 * Ny], vertices[i][2 * Ny]
                    });
            cell.index = n_cells;
            grid += cell;
            ++n_cells;
        }
    }

    grid.setup_adjacency();

    return grid;
}

} // namespace zephyr::geom::generator