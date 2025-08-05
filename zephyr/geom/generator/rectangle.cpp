#include <iostream>
#include <algorithm>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/rectangle.h>

#include <zephyr/geom/primitives/quad.h>
#include <zephyr/geom/primitives/polygon.h>

#include <zephyr/mesh/euler/amr_cells.h>

#include <zephyr/utils/json.h>

namespace zephyr::geom::generator {

using namespace zephyr::mesh;

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

void Rectangle::initialize(AmrCells& cells) {
    check_size();

    if (m_axial && m_ymin == 0.0) {
        m_bounds.bottom = Boundary::WALL;
    }

    if (m_voronoi) {
        return initialize_voronoi(cells);
    } else {
        return initialize_classic(cells);
    }
}

void Rectangle::initialize_classic(AmrCells& cells)  {
    bool x_period = periodic_along_x();
    bool y_period = periodic_along_y();

    double hx = (m_xmax - m_xmin) / m_nx;
    double hy = (m_ymax - m_ymin) / m_ny;

    auto get_index = [=](index_t i, index_t j) -> index_t {
        return i * m_ny + j;
    };

    auto get_index_pair = [=](index_t n) -> std::array<index_t, 2> {
        return {n / m_ny, n % m_ny};
    };

    auto get_vertex = [=](index_t i, index_t j) -> Vector3d {
        return {
                m_xmin + ((m_xmax - m_xmin) * i) / m_nx,
                m_ymin + ((m_ymax - m_ymin) * j) / m_ny,
                0.0
        };
    };

    auto neib_index = [=](index_t i, index_t j, Side2D side) -> index_t {
        if (side == Side2D::LEFT) {
            return i == 0 && !x_period ?  get_index(i, j) : get_index((i - 1 + m_nx) % m_nx, j);
        }
        else if (side == Side2D::RIGHT) {
            return i == m_nx - 1 && !x_period ? get_index(i, j) : get_index((i + 1 + m_nx) % m_nx, j);
        }
        else if (side == Side2D::BOTTOM) {
            return j == 0 && !y_period ? get_index(i, j) : get_index(i, (j - 1 + m_ny) % m_ny);
        }
        else if (side == Side2D::TOP) {
            return j == m_ny - 1 && !y_period ? get_index(i, j): get_index(i, (j + 1 + m_ny) % m_ny);
        }
        else {
            throw std::runtime_error("Strange side #142");
        }
    };

    cells.set_dimension(2);
    cells.set_adaptive(true);
    cells.set_linear(true);
    cells.set_axial(m_axial);

    cells.resize_amr(m_size);

    int n_faces = 8;
    int n_nodes = 9;

    for (index_t ic = 0; ic < m_size; ++ic) {
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

        cells.faces.boundary[iface + Side2D::L] = i > 0 ? Boundary::ORDINARY : m_bounds.left;
        cells.faces.boundary[iface + Side2D::R] = i < m_nx - 1 ? Boundary::ORDINARY : m_bounds.right;
        cells.faces.boundary[iface + Side2D::B] = j > 0 ? Boundary::ORDINARY : m_bounds.bottom;
        cells.faces.boundary[iface + Side2D::T] = j < m_ny - 1 ? Boundary::ORDINARY : m_bounds.top;

        for (auto side: Side2D::items()) {
            cells.faces.adjacent.rank[iface + side] = 0;
            cells.faces.adjacent.index[iface + side] = neib_index(i, j, side);
            cells.faces.adjacent.alien[iface + side] = -1;
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

        cells.faces.vertices[iface + Side2D::L] = Side2D::L.sf();
        cells.faces.vertices[iface + Side2D::R] = Side2D::R.sf();
        cells.faces.vertices[iface + Side2D::B] = Side2D::B.sf();
        cells.faces.vertices[iface + Side2D::T] = Side2D::T.sf();

        for (index_t jn = 0; jn < n_nodes; ++jn) {
            cells.verts[ic * n_nodes + jn] = quad[jn];
        }

        if (m_axial) {
            cells.volume_alt[ic] = hx * hy * quad.vs<0, 0>().y();
            cells.faces.area_alt[iface + Side2D::L] = hy * quad.vs<-1, 0>().y();
            cells.faces.area_alt[iface + Side2D::R] = hy * quad.vs<+1, 0>().y();
            cells.faces.area_alt[iface + Side2D::B] = hx * quad.vs< 0,-1>().y();
            cells.faces.area_alt[iface + Side2D::T] = hx * quad.vs< 0,+1>().y();
        }
    }
}

void Rectangle::initialize_voronoi(AmrCells& cells) {
    Grid grid = make();
    m_size = grid.n_cells();

    cells.set_dimension(2);
    cells.set_adaptive(false);
    cells.set_linear(true);
    cells.set_axial(m_axial);

    // Зарезервировать ячейки, добавляем push_back'ом
    cells.reserve(m_size, 6 * m_size, 6 * m_size);

    for (index_t ic = 0; ic < m_size; ++ic) {
        GCell cell = grid.cell(ic);

        int n_nodes = cell.n_nodes();

        // Уникальные узлы
        std::vector<GNode::Ptr> nodes(n_nodes, nullptr);

        Polygon poly(n_nodes);
        for (int i = 0; i < n_nodes; ++i) {
            nodes[i] = grid.node(cell.node(i).index);
            poly.set(i, nodes[i]->v);
        }

        cells.push_back(poly);

        // Смежность
        for (int i = 0; i < poly.size(); ++i) {
            index_t iface = cells.face_begin[ic] + i;

            GNode::Ref v1 = nodes[i];
            GNode::Ref v2 = nodes[(i + 1) % n_nodes];

            cells.faces.boundary[iface] = cell.boundary({v1, v2});
            cells.faces.adjacent.rank[iface]  = 0;
            cells.faces.adjacent.index[iface] = cell.adjacent({v1, v2});
            cells.faces.adjacent.alien[iface] = -1;
            cells.faces.adjacent.basic[iface] = ic;
        }
    }
}

} // namespace zephyr::geom::generator