#include <iostream>
#include <algorithm>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/primitives/quad.h>
#include <zephyr/utils/json.h>
#include <zephyr/geom/generator/strip.h>
#include <zephyr/mesh/euler/amr_cells.h>

namespace zephyr::geom::generator {

using namespace zephyr::mesh;

Strip::Strip(const Json& config)
        : Generator("strip"),
          m_xmin(0.0), m_xmax(1.0), m_nx(0) {

    if (!config["geometry"]) {
        throw std::runtime_error("Strip config doesn't contain key 'geometry'");
    }

    m_xmin = config["geometry"]["x_min"].as<double>();
    m_xmax = config["geometry"]["x_max"].as<double>();
    
    if (!config["bounds"]) {
        throw std::runtime_error("Strip config doesn't contain key 'bounds'");
    }
    m_bounds.left   = boundary_from_string(config["bounds"]["left"].as<std::string>());
    m_bounds.right  = boundary_from_string(config["bounds"]["right"].as<std::string>());

    if (!config["size"] || !config["size"].is_number()) {
        throw std::runtime_error("Strip config doesn't contain key 'size'");
    }

    set_size(config["cells"].as<int>());
}

Strip::Strip(double xmin, double xmax, Type type) :
        Generator("strip"),
        m_type(type),
        m_nx(0),
        m_xmin(xmin),
        m_xmax(xmax),
        m_bounds() {
    check_params();
}

int Strip::size() const {
    return m_nx;
}

Box Strip::bbox() const {
    Vector3d vmin(x_min(), y_min(), 0.0);
    Vector3d vmax(x_max(), y_max(), 0.0);

    return {vmin, vmax};
}

void Strip::set_nx(int nx) {
    set_size(nx);
}

void Strip::set_size(int N) {
    if (N < 1) {
        std::cerr << "Strip Error: N < 1\n";
        throw std::runtime_error("Strip Error: N < 1");
    }
    if (N > 1000000000) {
        std::cerr << "Attempt to create mesh with more than 1 billion cells\n";
        throw std::runtime_error("Attempt to create mesh with more than 1 billion cells");
    }

    m_nx = N;
}

void Strip::set_boundaries(Boundaries bounds) {
    m_bounds = bounds;
    if (periodic_along_x()) {
        m_bounds.left = m_bounds.right = Boundary::PERIODIC;
    }
}

double Strip::x_min() const {
    return m_xmin;
}

double Strip::x_max() const {
    return m_xmax;
}

double Strip::y_min() const {
    return -0.5 * aspect * (x_max() - x_min());
}

double Strip::y_max() const {
    return +0.5 * aspect * (x_max() - x_min());
}

int Strip::nx() const {
    return m_nx;
}

bool Strip::periodic_along_x() const {
    return m_bounds.left == Boundary::PERIODIC || m_bounds.right == Boundary::PERIODIC;
}

void Strip::check_params() const {
    if (m_xmin >= m_xmax) {
        std::cerr << "Strip Error: x_min >= x_max\n";
        throw std::runtime_error("Strip Error: x_min >= x_max");
    }
}

std::vector<double> nodes_uniform(double xmin, double xmax, int size) {
    std::vector<double> nodes(size + 1);
    for (int i = 0; i <= size; ++i) {
        nodes[i] = xmin + (xmax - xmin) * i / size;
    }
    return nodes;
}

std::vector<double> nodes_random(double xmin, double xmax, int size) {
    std::vector<double> nodes(size + 1);
    for (int i = 0; i <= size; ++i) {
        nodes[i] = xmin + (xmax - xmin) * rand() / double(RAND_MAX);
    }

    std::sort(nodes.begin(), nodes.end());
    nodes[0] = xmin;
    nodes[size] = xmax;

    for (int i = 1; i < size; ++i) {
        nodes[i] = (nodes[i - 1] + 2 * nodes[i] + nodes[i + 1]) / 4.0;
    }
    return nodes;
}

std::vector<double> get_nodes(Strip::Type type, double xmin, double xmax, int size) {
    if (type == Strip::Type::UNIFORM) {
        return nodes_uniform(xmin, xmax, size);
    }
    else {
        return nodes_random(xmin, xmax, size);
    }
}

Grid Strip::make() {
    check_size();

    auto nodes1D = get_nodes(m_type, m_xmin, m_xmax, m_nx);

    double y1 = y_min();
    double y2 = y_max();

    std::vector<std::vector<GNode::Ptr>> nodes(
            2, std::vector<GNode::Ptr>(m_nx + 1, nullptr));

    Grid grid;

    grid.reserve_nodes(m_nx + 1);
    int n_nodes = 0;
    for (int i = 0; i <= m_nx; ++i) {
        nodes[0][i] = GNode::create(nodes1D[i], y1);
        nodes[0][i]->index = n_nodes;
        grid += nodes[0][i];
        ++n_nodes;

        nodes[1][i] = GNode::create(nodes1D[i], y2);
        nodes[1][i]->index = n_nodes;
        grid += nodes[1][i];
        ++n_nodes;
    }
    nodes[0][0]->add_boundary(m_bounds.left);
    nodes[1][0]->add_boundary(m_bounds.left);
    nodes[0][m_nx]->add_boundary(m_bounds.right);
    nodes[1][m_nx]->add_boundary(m_bounds.right);

    grid.reserve_cells(m_nx);
    for (int i = 0; i < m_nx; ++i) {
        GCell cell = GCell::quad(
                {
                        nodes[0][i], nodes[0][i + 1],
                        nodes[1][i + 1], nodes[1][i]
                });
        cell.index = i;
        grid += cell;
    }

    grid.setup_adjacency();

    return grid;
}

void Strip::initialize(AmrCells& cells) {
    bool x_period = periodic_along_x();

    double m_ymin = y_min();
    double m_ymax = y_max();

    index_t m_ny = 1;

    double hx = (m_xmax - m_xmin) / m_nx;
    double hy = (m_ymax - m_ymin) / m_ny;

    auto get_vertex = [=](index_t i, index_t j) -> Vector3d {
        return {
                m_xmin + ((m_xmax - m_xmin) * i) / m_nx,
                m_ymin + ((m_ymax - m_ymin) * j) / m_ny,
                0.0
        };
    };

    auto neib_index = [=](index_t i, Side2D side) -> index_t {
        if (side == Side2D::LEFT) {
            return i == 0 && !x_period ?  i : (i - 1 + m_nx) % m_nx;
        }
        else if (side == Side2D::RIGHT) {
            return i == m_nx - 1 && !x_period ? i : (i + 1 + m_nx) % m_nx;
        }
        else {
            throw std::runtime_error("Strange side #153");
        }
    };

    cells.set_dimension(2);
    cells.set_adaptive(true);
    cells.set_linear(true);
    cells.set_axial(false);

    cells.resize(m_nx);

    int n_faces = 8;
    int n_nodes = 9;

    for (index_t ic = 0; ic < m_nx; ++ic) {
        cells.next[ic] = ic;
        cells.rank[ic] = 0;
        cells.index[ic] = ic;

        cells.flag[ic] = 0;
        cells.level[ic] = 0;
        cells.b_idx[ic] = ic;
        cells.z_idx[ic] = 0;

        SqQuad quad(
                get_vertex(ic, 0),
                get_vertex(ic + 1, 0),
                get_vertex(ic, 1),
                get_vertex(ic + 1, 1));

        cells.center[ic] = quad.vs<0, 0>();
        cells.volume[ic] = hx * hy;
        cells.volume_alt[ic] = NAN;
        cells.face_begin[ic] = ic * n_faces;
        cells.node_begin[ic] = ic * n_nodes;
        cells.face_begin[ic + 1] = cells.face_begin[ic] + n_faces;
        cells.node_begin[ic + 1] = cells.face_begin[ic] + n_nodes;

        // INIT FACES
        for (index_t iface: cells.faces_range(ic)) {
            cells.faces.set_undefined(iface);
            cells.faces.area_alt[iface] = NAN;
        }

        index_t iface = ic * n_faces;

        cells.faces.boundary[iface + Side2D::L] = ic > 0 ? Boundary::ORDINARY : m_bounds.left;
        cells.faces.boundary[iface + Side2D::R] = ic < m_nx - 1 ? Boundary::ORDINARY : m_bounds.right;

        for (auto side: {Side2D::LEFT, Side2D::RIGHT}) {
            cells.faces.adjacent.rank[iface + side] = 0;
            cells.faces.adjacent.index[iface + side] = neib_index(ic, side);
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
    }
}

} // namespace zephyr::geom::generator