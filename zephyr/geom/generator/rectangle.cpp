#include <iostream>
#include <algorithm>
#include <format>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/primitives/quad.h>
#include <zephyr/utils/json.h>
#include <zephyr/mesh/side.h>

namespace zephyr::geom::generator {

using mesh::Side2D;

Rectangle::Rectangle(const Json& config)
    : Generator("rectangle"),
      m_xmin(0.0), m_xmax(1.0),
      m_ymin(0.0), m_ymax(1.0) {

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
            throw std::runtime_error(std::format("Rectangle config strange sizes: {}, {}, "
                                                 "setup key size.nx, size.ny or both", nx, ny));
        }
    }
}

Rectangle::Rectangle()
    : Rectangle(0.0, 1.0, 0.0, 1.0, false) {
    set_boundaries(Boundaries{});
}

Rectangle::Rectangle(double xmin, double xmax, double ymin, double ymax, bool voronoi) :
        Generator("rectangle"),
        m_xmin(xmin), m_xmax(xmax),
        m_ymin(ymin), m_ymax(ymax),
        m_voronoi(voronoi) {
    check_params();
}

Box Rectangle::bbox() const {
    Vector3d vmin(m_xmin, m_ymin, 0.0);
    Vector3d vmax(m_xmax, m_ymax, 0.0);

    return {vmin, vmax};
}

void Rectangle::set_axial(bool axial) {
    m_axial = axial;

    if (m_axial && m_ymin == 0.0) {
        m_bounds.bottom = Boundary::WALL;
    }
}

void Rectangle::set_nx(int nx) {
    if (nx < 1) {
        throw std::runtime_error("Rectangle: Nx < 1");
    }
    if (!m_voronoi) {
        m_nx = nx;
        m_ny = std::max(int(std::round(m_nx * (m_ymax - m_ymin) / (m_xmax - m_xmin))), 1);
    }
    else {
        m_nx = std::max(int(std::round(nx * std::pow(3.0, 0.25)) / std::sqrt(2.0)), 1);
        m_ny = std::max(int(std::round(m_nx * (m_ymax - m_ymin) / (m_xmax - m_xmin) / std::sqrt(3.0))), 1);
    }
    compute_size();
}

void Rectangle::set_ny(int ny) {
    if (ny < 1) {
        throw std::runtime_error("Rectangle: Ny < 1");
    }
    if (!m_voronoi) {
        m_ny = ny;
        m_nx = std::max(int(std::round(ny * (m_xmax - m_xmin) / (m_ymax - m_ymin))), 1);
    }
    else {
        m_ny = std::max(int(std::round(ny / std::pow(3.0, 0.25)) / std::sqrt(2.0)), 1);
        m_nx = std::max(int(std::round(m_ny * (m_xmax - m_xmin) / (m_ymax - m_ymin) * std::sqrt(3.0))), 1);
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
        m_nx = std::max(int(std::round(nx * std::pow(3.0, 0.25)) / std::sqrt(2.0)), 1);
        m_ny = std::max(int(std::round(ny / std::pow(3.0, 0.25)) / std::sqrt(2.0)), 1);
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
        m_nx = std::max(int(std::round((m_xmax - m_xmin) / d)), 1);
        m_ny = std::max(int(std::round((m_ymax - m_ymin) / d)), 1);
    }
    else {
        double a = std::sqrt(2 * (m_xmax - m_xmin) * (m_ymax - m_ymin) / (std::sqrt(3.0) * N));
        double h = std::sqrt(3.0) * a;

        m_nx = std::max(int(std::round((m_xmax - m_xmin) / a)), 1);
        m_ny = std::max(int(std::round((m_ymax - m_ymin) / h)), 1);
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
    if (m_axial && m_ymin == 0.0) {
        m_bounds.bottom = Boundary::WALL;
    }
}

void Rectangle::set_adaptive(bool adaptive) {
    if (!m_voronoi) {
        m_adaptive = adaptive;
    }
    else {
        m_adaptive = false;
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
        throw std::runtime_error("Rectangle::check_params: x_min >= x_max");
    }
    if (m_ymin >= m_ymax) {
        throw std::runtime_error("Rectangle::check_params: y_min >= y_max");
    }
}

void Rectangle::compute_size() {
    if (!m_voronoi) {
        m_size = m_nx * m_ny;
    } else {
        m_size = 2 * m_nx * m_ny;
    }
    if (m_size > max_grid_size) {
        throw std::runtime_error(std::format("Generator::check_size: attempt to create mesh "
                                             "that contains more than {} elements", max_grid_size));
    }
}

Grid Rectangle::make() const {
    if (m_voronoi) {
        return create_voronoi();
    }
    if (m_adaptive) {
        return create_classic_amr();
    }
    return create_classic();
}

Grid Rectangle::create_classic() const {
    check_size(m_size);

    double dx = (m_xmax - m_xmin) / m_nx;
    double dy = (m_ymax - m_ymin) / m_ny;

    Grid grid;
    grid.reserve_nodes((m_nx + 1) * (m_ny + 1));
    grid.reserve_cells(m_nx * m_ny);

    std::vector nodes(m_nx + 1, std::vector<NodeInput::Ptr>(m_ny + 1));
    for (int i = 0; i <= m_nx; ++i) {
        for (int j = 0; j <= m_ny; ++j) {
            double x = m_xmin + i * dx;
            double y = m_ymin + j * dy;
            nodes[i][j] = NodeInput::create({x, y, 0.0});

            grid.add_node(nodes[i][j]);
        }
    }

    // Ячейки как полигоны, обход граней и вершин против часовой
    // стрелки, начиная с нижней левой вершины (нижней грани)
    std::array<Boundary, 4> bc{};
    for (int i = 0; i < m_nx; ++i) {
        bc[3] = i == 0      ? m_bounds.left   : Boundary::INNER;
        bc[1] = i == m_nx-1 ? m_bounds.right  : Boundary::INNER;
        for (int j = 0; j < m_ny; ++j) {
            bc[0] = j == 0      ? m_bounds.bottom : Boundary::INNER;
            bc[2] = j == m_ny-1 ? m_bounds.top    : Boundary::INNER;

            grid.add_quad({
                nodes[i][j],
                nodes[i + 1][j],
                nodes[i + 1][j + 1],
                nodes[i][j + 1]},
                bc);
        }
    }
    return grid;
}


Grid Rectangle::create_classic_amr() const {
    check_size(m_size);

    double dx = (m_xmax - m_xmin) / m_nx;
    double dy = (m_ymax - m_ymin) / m_ny;

    Grid grid;
    grid.reserve_nodes((m_nx + 1) * (m_ny + 1));
    grid.reserve_cells(m_nx * m_ny);

    std::vector nodes(m_nx + 1, std::vector<NodeInput::Ptr>(m_ny + 1));
    for (int i = 0; i <= m_nx; ++i) {
        for (int j = 0; j <= m_ny; ++j) {
            double x = m_xmin + i * dx;
            double y = m_ymin + j * dy;
            nodes[i][j] = NodeInput::create({x, y, 0.0});
            grid.add_node(nodes[i][j]);
        }
    }

    // AMR-ячейки, обход граней как в Side2D, вершины в Z-порядке
    std::vector<Boundary> bc(4);
    for (int i = 0; i < m_nx; ++i) {
        bc[Side2D::L] = i == 0      ? m_bounds.left   : Boundary::INNER;
        bc[Side2D::R] = i == m_nx-1 ? m_bounds.right  : Boundary::INNER;
        for (int j = 0; j < m_ny; ++j) {
            bc[Side2D::B] = j == 0      ? m_bounds.bottom : Boundary::INNER;
            bc[Side2D::T] = j == m_ny-1 ? m_bounds.top    : Boundary::INNER;
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
    check_size(m_size);

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
    std::vector vertices(Nx + 2, std::vector<NodeInput::Ptr>(2 * Ny + 1, nullptr));

    for (size_t j = 0; j <= 2 * Ny; ++j) {
        double y = y_shift + h * j;

        // Часть вершин на левой границе пропускаем
        if (j % 2 == 0) {
            vertices[0][j] = NodeInput::create({m_xmin, y, 0.0});
            vertices[0][j]->bc = m_bounds.left;
        }

        for (size_t i = 1; i <= Nx; ++i) {
            double x = x_shift + double(3 * i - (i + j) % 2) * 0.5 * D;
            vertices[i][j] = NodeInput::create({x, y, 0.0});

            if (j == 0) {
                vertices[i][j]->bc = m_bounds.bottom;
            } else if (j == 2 * Ny) {
                vertices[i][j]->bc = m_bounds.top;
            }
        }

        // Часть вершин на правой границе пропускаем
        if (j == 0 || j == 2 * Ny || j % 2 == Nx % 2) {
            vertices[Nx + 1][j] = NodeInput::create({m_xmax, y, 0.0});
            vertices[Nx + 1][j]->bc = m_bounds.right;
        }
    }

    Grid grid;

    // Переносим вершины в один массив
    grid.reserve_nodes((Nx + 2)*(2 * Ny + 1));
    for (size_t i = 0; i <= Nx + 1; ++i) {
        for (size_t j = 0; j <= 2 * Ny; ++j) {
            if (!vertices[i][j]->pos.hasNaN()) {
                grid.add_node(vertices[i][j]);
            }
        }
    }

    using VList = std::vector<NodeInput::Ptr>;

    auto erase_nans = [](VList& vlist) {
        const auto to_remove = std::ranges::remove_if(vlist,
            [](NodeInput::Ref v) -> bool { return !v; }).begin();
        vlist.erase(to_remove, vlist.end());
    };

    // Создаем массив ячеек
    grid.reserve_nodes(m_size);
    for (size_t i = 0; i <= Nx; ++i) {
        if (i % 2 == 0) {
            for (size_t j = 0; j < 2 * Ny - 1; j += 2) {
                VList vlist = {
                    vertices[i][j], vertices[i + 1][j],
                    vertices[i + 1][j + 1], vertices[i + 1][j + 2],
                    vertices[i][j + 2], vertices[i][j + 1]
                };
                erase_nans(vlist);

                grid.add_cell(CellType::POLYGON, vlist);
            }
        } else {
            // j == 0
            grid.add_cell(CellType::QUAD, {
                vertices[i][0], vertices[i + 1][0],
                vertices[i + 1][1], vertices[i][1]
            });

            for (size_t j = 1; j < 2 * Ny - 2; j += 2) {
                VList vlist = {
                    vertices[i][j], vertices[i + 1][j],
                    vertices[i + 1][j + 1], vertices[i + 1][j + 2],
                    vertices[i][j + 2], vertices[i][j + 1]
                };
                erase_nans(vlist);
                grid.add_cell(CellType::POLYGON, vlist);
            }

            // j == last; i % 2 == 1
            grid.add_cell(CellType::QUAD, {
                vertices[i][2 * Ny - 1], vertices[i + 1][2 * Ny - 1],
                vertices[i + 1][2 * Ny], vertices[i][2 * Ny]
            });
        }
    }

    return grid;
}

} // namespace zephyr::geom::generator