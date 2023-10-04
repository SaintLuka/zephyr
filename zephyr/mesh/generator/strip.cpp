#include <iostream>
#include <algorithm>

#include <zephyr/utils/mpi.h>
#include <zephyr/geom/cell.h>

#include <zephyr/geom/box.h>
#include <zephyr/mesh/generator/strip.h>


namespace zephyr { namespace mesh { namespace generator {

using namespace zephyr::geom;
using zephyr::utils::mpi;

Strip::Strip(double xmin, double xmax, Nodes nodes) :
        Generator("strip"),
        m_nodes(nodes), m_xmin(xmin), m_xmax(xmax),
        m_left_flag(FaceFlag::WALL), m_right_flag(FaceFlag::WALL)
        {
    check_params();
}

Box Strip::bbox() const {
    Vector3d vmin(x_min(), y_min(), 0.0);
    Vector3d vmax(x_max(), y_max(), 0.0);

    return Box(vmin, vmax);
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

    m_size = N;
}

void Strip::set_boundary_flags(FaceFlag left, FaceFlag right) {
    m_left_flag = left;
    m_right_flag = right;
    if (periodic_along_x()) {
        m_left_flag = m_right_flag = FaceFlag::PERIODIC;
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
    return m_size;
}

bool Strip::periodic_along_x() const {
    return m_left_flag == FaceFlag::PERIODIC || m_right_flag == FaceFlag::PERIODIC;
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

std::vector<double> nodes_sine(double xmin, double xmax, int size) {
    std::vector<double> nodes(size + 1);

    for (int i = 0; i <= size; ++i) {
        nodes[i] = xmin + (xmax - xmin) * i / size;
    }

    return nodes;
}

std::vector<double> nodes_exp(double xmin, double xmax, int size) {
    std::vector<double> nodes(size + 1);

    for (int i = 0; i <= size; ++i) {
        nodes[i] = xmin + (xmax - xmin) * i / size;
    }

    return nodes;
}

std::vector<double> get_nodes(Strip::Nodes type, double xmin, double xmax, int size) {
    switch (type) {
        case Strip::Nodes::SINE:
            return nodes_sine(xmin, xmax, size);
        case Strip::Nodes::EXP:
            return nodes_exp(xmin, xmax, size);
        case Strip::Nodes::RANDOM:
            return nodes_random(xmin, xmax, size);
        default:
            return nodes_uniform(xmin, xmax, size);
    }
}

void Strip::initialize(Storage &cells) {
    cells.resize(m_size);

    auto nodes = get_nodes(m_nodes, m_xmin, m_xmax, m_size);

    double y1 = y_min();
    double y2 = y_max();

    for (int i = 0; i < m_size; ++i) {
        // Границы ячейки
        double x1 = nodes[i];
        double x2 = nodes[i+1];

        ShortList2D verts = {
                Vector3d(x1, y1, 0.0), Vector3d(x2, y1, 0.0),
                Vector3d(x1, y2, 0.0), Vector3d(x2, y2, 0.0)
        };

        Cell cell(verts);

        // Базовая информация
        cell.rank  = mpi::rank();
        cell.index = i;

        // Данные AMR
        cell.b_idx = i;
        cell.z_idx = 0;
        cell.next  = 0;
        cell.level = 0;
        cell.flag  = 0;

        // Флаги граничных условий
        auto ordinary = FaceFlag::ORDINARY;

        cell.faces[Side::L].boundary = i > 0 ? ordinary : m_left_flag;
        cell.faces[Side::R].boundary = i < m_size - 1 ? ordinary : m_right_flag;
        cell.faces[Side::B].boundary = FaceFlag::UNDEFINED;
        cell.faces[Side::T].boundary = FaceFlag::UNDEFINED;

        cell.faces[Side::L].adjacent.rank = 0;
        cell.faces[Side::R].adjacent.rank = 0;
        cell.faces[Side::B].adjacent.rank = 0;
        cell.faces[Side::T].adjacent.rank = 0;

        cell.faces[Side::L].adjacent.ghost = -1;
        cell.faces[Side::R].adjacent.ghost = -1;
        cell.faces[Side::B].adjacent.ghost = -1;
        cell.faces[Side::T].adjacent.ghost = -1;

        cell.faces[Side::L].adjacent.index = std::max(0, i - 1);
        cell.faces[Side::R].adjacent.index = std::min(i + 1, m_size - 1);
        cell.faces[Side::B].adjacent.index = i;
        cell.faces[Side::T].adjacent.index = i;

        cells[i].geom() = cell;
    }
}

} // generator
} // mesh
} // zephyr