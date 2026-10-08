#include <iostream>
#include <algorithm>
#include <random>

#include <zephyr/geom/side.h>
#include <zephyr/geom/boundary.h>
#include <zephyr/geom/indexing.h>
#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/primitives/quad.h>
#include <zephyr/geom/generator/strip.h>
#include <zephyr/utils/json.h>
#include <zephyr/mesh/raw/raw_cells.h>

namespace zephyr::geom::generator {

using namespace zephyr::mesh;

Strip::Strip(const Json& config)
    : Generator("strip"), x_min_(0.0), x_max_(1.0) {

    if (!config["geometry"]) {
        throw std::runtime_error("Strip config doesn't contain key 'geometry'");
    }

    x_min_ = config["geometry"]["x_min"].as<double>();
    x_max_ = config["geometry"]["x_max"].as<double>();
    
    if (!config["bounds"]) {
        throw std::runtime_error("Strip config doesn't contain key 'bounds'");
    }
    bounds_.left  = boundary_from_string(config["bounds"]["left"].as<std::string>());
    bounds_.right = boundary_from_string(config["bounds"]["right"].as<std::string>());

    if (!config["size"] || !config["size"].is_number()) {
        throw std::runtime_error("Strip config doesn't contain key 'size'");
    }

    set_size(config["cells"].as<int>());
}

Strip::Strip(double x_min, double x_max, Type type) :
        Generator("strip"),
        m_type(type),
        x_min_(x_min),
        x_max_(x_max) {
    check_params();
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
        throw std::runtime_error("Strip::set_size: N < 1");
    }
    if (N > 200'000'000) {
        std::cerr << "Attempt to create mesh with more than 1 billion cells\n";
        throw std::runtime_error("Attempt to create mesh with more than 1 billion cells");
    }

    nx_ = N;
}

void Strip::set_boundaries(Boundaries bounds) {
    bounds_ = bounds;
    if (periodic_along_x()) {
        bounds_.left = bounds_.right = Boundary::PERIODIC;
    }
}

double Strip::x_min() const {
    return x_min_;
}

double Strip::x_max() const {
    return x_max_;
}

double Strip::y_min() const {
    return -0.5 * aspect * (x_max() - x_min());
}

double Strip::y_max() const {
    return +0.5 * aspect * (x_max() - x_min());
}

int Strip::nx() const {
    return nx_;
}

bool Strip::periodic_along_x() const {
    return bounds_.left == Boundary::PERIODIC || bounds_.right == Boundary::PERIODIC;
}

void Strip::check_params() const {
    if (x_min_ >= x_max_) {
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
    static std::default_random_engine gen;
    std::uniform_real_distribution uniform(xmin, xmax);

    std::vector<double> nodes(size + 1);
    for (int i = 0; i <= size; ++i) {
        nodes[i] = uniform(gen);
    }

    std::ranges::sort(nodes);
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

Grid Strip::make() const {
    check_size(nx_);

    auto nodes1D = get_nodes(m_type, x_min_, x_max_, nx_);

    double y1 = y_min();
    double y2 = y_max();

    std::vector nodes(2, std::vector<GNode::Ptr>(nx_ + 1, nullptr));

    Grid grid;

    grid.reserve_nodes(nx_ + 1);
    for (int i = 0; i <= nx_; ++i) {
        nodes[0][i] = GNode::create({nodes1D[i], y1, 0.0});
        grid.add_node(nodes[0][i]);

        nodes[1][i] = GNode::create({nodes1D[i], y2, 0.0});
        grid.add_node(nodes[1][i]);
    }
    nodes[0][0]->bc = bounds_.left;
    nodes[1][0]->bc = bounds_.left;
    nodes[0][nx_]->bc = bounds_.right;
    nodes[1][nx_]->bc = bounds_.right;

    grid.reserve_cells(nx_);
    for (int i = 0; i < nx_; ++i) {
        grid.add_cell(
            CellType::QUAD, {
                nodes[0][i], nodes[0][i + 1],
                nodes[1][i + 1], nodes[1][i]
            });
    }

    return grid;
}

RawCells Strip::make_cells(bool unique_nodes) const {
    bool x_period = periodic_along_x();

    double y_min_ = y_min();
    double y_max_ = y_max();

    index_t ny_ = 1;

    double hx = (x_max_ - x_min_) / nx_;
    double hy = (y_max_ - y_min_) / ny_;

    auto get_vertex = [=, this](index_t i, index_t j) -> Vector3d {
        return {
                x_min_ + ((x_max_ - x_min_) * i) / nx_,
                y_min_ + ((y_max_ - y_min_) * j) / ny_,
                0.0
        };
    };

    auto neib_index = [=, this](index_t i, Side2D side) -> index_t {
        if (side == Side2D::LEFT) {
            return i == 0 && !x_period ?  i : (i - 1 + nx_) % nx_;
        }
        if (side == Side2D::RIGHT) {
            return i == nx_ - 1 && !x_period ? i : (i + 1 + nx_) % nx_;
        }
        throw std::runtime_error("Strange side #153");
    };

    RawCells cells({
        .dim = 2,
        .adaptive = true,
        .linear = true,
        .axial = false,
        .nodes = unique_nodes
    });

    cells.resize_amr(nx_);

    static constexpr int n_faces = 8;
    static constexpr int n_nodes = 9;

    for (index_t ic = 0; ic < nx_; ++ic) {
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
        cells.verts.offsets[ic] = ic * n_nodes;
        cells.verts.offsets[ic + 1] = cells.verts.offsets[ic] + n_nodes;

        // INIT FACES
        cells.faces.offsets[ic] = ic * n_faces;
        cells.faces.offsets[ic + 1] = cells.faces.offsets[ic] + n_faces;

        for (index_t iface: cells.faces.range(ic)) {
            cells.faces.set_undefined(iface);
            cells.faces.area_alt[iface] = NAN;
        }

        index_t iface = ic * n_faces;

        cells.faces.boundary[iface + Side2D::L] = ic > 0 ? Boundary::INNER : bounds_.left;
        cells.faces.boundary[iface + Side2D::R] = ic < nx_ - 1 ? Boundary::INNER : bounds_.right;

        for (auto side: {Side2D::LEFT, Side2D::RIGHT}) {
            cells.faces.adjacent.rank[iface + side] = 0;
            cells.faces.adjacent.index[iface + side] = neib_index(ic, side);
            cells.faces.adjacent.ghost[iface + side] = -1;
            cells.faces.adjacent.basic[iface + side] = ic;
            cells.faces.vertices[iface + side].fill(-1);
        }

        for (auto side: {Side2D::BOTTOM, Side2D::TOP}) {
            cells.faces.adjacent.rank[iface + side] = 0;
            cells.faces.adjacent.index[iface + side] = ic;
            cells.faces.adjacent.ghost[iface + side] = -1;
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

        cells.faces.vertices[iface + Side2D::L] = indexing::amr::sf(Side2D::L);
        cells.faces.vertices[iface + Side2D::R] = indexing::amr::sf(Side2D::R);
        cells.faces.vertices[iface + Side2D::B] = indexing::amr::sf(Side2D::B);
        cells.faces.vertices[iface + Side2D::T] = indexing::amr::sf(Side2D::T);

        for (index_t jn = 0; jn < n_nodes; ++jn) {
            cells.verts[ic * n_nodes + jn] = quad[jn];
        }
    }
    return cells;
}

} // namespace zephyr::geom::generator