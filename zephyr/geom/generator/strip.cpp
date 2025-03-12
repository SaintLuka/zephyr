#include <iostream>
#include <algorithm>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/utils/json.h>
#include <zephyr/geom/generator/strip.h>

namespace zephyr::geom::generator {

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

} // namespace zephyr::geom::generator