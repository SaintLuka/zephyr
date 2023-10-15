#include <cassert>
#include <zephyr/utils/mpi.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/primitives/amr_cell.h>

namespace zephyr::geom {

using namespace zephyr::utils;


GNode::GNode(const Vector3d &_v)
    : index(-1), v(_v) { }

GNode::GNode(double x, double y, double z)
    : index(-1), v({x, y, z}) { }

GNode::Ptr GNode::create(const Vector3d &v) {
    return std::make_shared<GNode>(v);
}

GNode::Ptr GNode::create(double x, double y, double z) {
    return std::make_shared<GNode>(x, y, z);
}

void GNode::add_boundary(Boundary flag) {
    m_bounds.insert(flag);
}

Boundary GNode::face_boundary(const std::vector<GNode::Ptr>& vs) {
    if (vs.empty()) {
        return Boundary::UNDEFINED;
    }

    std::set<Boundary> intersection = vs[0]->m_bounds;
    for (size_t i = 1; i < vs.size(); ++i) {
        std::set<Boundary> inter;

        std::set_intersection(
                intersection.begin(), intersection.end(),
                vs[i]->m_bounds.begin(), vs[i]->m_bounds.end(),
                std::inserter(inter, inter.begin()));

        intersection = inter;
    }

    if (intersection.empty()) {
        return Boundary::UNDEFINED;
    }
    else {
        return *intersection.begin();
    }
}

std::vector<int> GNode::shared_cells(const std::vector<GNode::Ptr> &vs) {
    if (vs.empty()) {
        return {};
    }

    assert(!vs[0]->m_neibs.empty());

    std::set<int> intersection = vs[0]->m_neibs;
    for (size_t i = 1; i < vs.size(); ++i) {
        assert(!vs[i]->m_neibs.empty());

        std::set<int> inter;

        std::set_intersection(
                intersection.begin(), intersection.end(),
                vs[i]->m_neibs.begin(), vs[i]->m_neibs.end(),
                std::inserter(inter, inter.begin()));

        intersection = inter;
    }

    if (intersection.empty()) {
        return {};
    }
    else {
        std::vector<int> res;
        res.reserve(intersection.size());
        for (auto idx: intersection) {
            res.push_back(idx);
        }
        return res;
    }
}

void GNode::add_neib_cell(int idx) {
    m_neibs.insert(idx);
}

GCell::GCell()
    : index(-1) { }

GCell::GCell(CellType type)
    : index(-1), m_type(type) { }

GCell GCell::triangle(const std::array<GNode::Ptr, 3>& nodes) {
    GCell cell(CellType::TRIANGLE);
    cell.m_nodes = {nodes[0], nodes[1], nodes[2]};
    cell.m_faces = {
            {0, 1},
            {1, 2},
            {2, 0}
    };
    return cell;
}

GCell GCell::quad(const std::array<GNode::Ptr, 4>& nodes) {
    GCell cell(CellType::QUAD);

    // Проверка нумерации вершин в m_nodes
    Vector3d a = nodes[0]->v - nodes[1]->v;
    Vector3d b = nodes[2]->v - nodes[1]->v;
    Vector3d c = nodes[2]->v - nodes[3]->v;
    Vector3d d = nodes[0]->v - nodes[3]->v;

    if (b.cross(a).dot(c.cross(d)) > 0.0) {
        cell.m_nodes = {nodes[0], nodes[1], nodes[2], nodes[3]};
    } else {
        cell.m_nodes = {nodes[0], nodes[1], nodes[3], nodes[2]};
    }
    cell.m_faces = {
            {0, 3},
            {1, 2},
            {0, 1},
            {3, 2}
    };
    cell.m_neibs = {-1, -1, -1, -1};
    return cell;
}

GCell GCell::polygon(const std::vector<GNode::Ptr>& nodes) {
    throw std::runtime_error("POLY");
}

GCell GCell::tetra(const std::array<GNode::Ptr, 4>& nodes) {
    throw std::runtime_error("TETRA");
}

GCell GCell::pyramid(const std::array<GNode::Ptr, 5>& nodes) {
    throw std::runtime_error("PYRAMID");
}

GCell GCell::wedge(const std::array<GNode::Ptr, 6>& nodes) {
    throw std::runtime_error("WEDGE");
}

GCell GCell::hexagedron(const std::array<GNode::Ptr, 6>& nodes) {
    throw std::runtime_error("HEX");
}

CellType GCell::type() const {
    return m_type;
}

GNode& GCell::node(int idx) {
    assert(idx < m_nodes.size());
    return *m_nodes[idx];
}

const GNode& GCell::node(int idx) const {
    assert(idx < m_nodes.size());
    return *m_nodes[idx];
}

int GCell::adjacent(int side) const {
    assert(side < m_neibs.size());
    return m_neibs[side];
}

void GCell::add_self_to_nodes() {
    assert(index >= 0);
    for (auto& node: m_nodes) {
        node->add_neib_cell(index);
    }
}

int GCell::adjacent(const std::vector<GNode::Ptr>& face_nodes) const {
    assert(index >= 0);

    auto neibs = GNode::shared_cells(face_nodes);

    if (neibs.empty()) {
        throw std::runtime_error("Grid::setup_adjacency error #1");
    }
    else if (neibs.size() == 1) {
        // Вероятно граничное условие
        if (neibs[0] == index) {
            // Все нормально граничное условие
            return index;
        }
        else {
            // Грань ссылается на непонятную ячейку
            throw std::runtime_error("Grid::setup_adjacency error #2");
        }
    }
    else if (neibs.size() == 2) {
        // Нормальные ситуации: по два соседа, при этом один из
        // соседей является самой ячейкой
        if (neibs[0] == index && neibs[1] != index) {
            return neibs[1];
        }
        else if (neibs[1] == index && neibs[0] != index) {
            return neibs[0];
        }
        else {
            // Два соседа, ни один не ссылается на саму ячеку
            throw std::runtime_error("Grid::setup_adjacency error #3");
        }
    }
    else {
        // Грань ссылается на множество ячеек
        throw std::runtime_error("Grid::setup_adjacency error #4");
    }
}

void GCell::find_neibs() {
    assert(index >= 0);
    for (int j = 0; j < n_faces(); ++j) {
        m_neibs[j] = adjacent(face_nodes(j));
    }
}

int GCell::n_nodes() const {
    return m_nodes.size();
}

int GCell::n_faces() const {
    return m_faces.size();
}

std::vector<GNode::Ptr> GCell::face_nodes(int idx) const {
    assert(idx < m_faces.size());
    std::vector<GNode::Ptr> nodes(m_faces[idx].size());
    for (int iv = 0; iv < int(m_faces[idx].size()); ++iv) {
        nodes[iv] = m_nodes[iv];
    }
    return nodes;
}

Boundary GCell::boundary(const std::vector<GNode::Ptr>& face_nodes) const {
    assert(index >= 0);

    int adj = adjacent(face_nodes);
    if (adj != index) {
        // Внутренняя грань
        return Boundary::ORDINARY;
    }
    else {
        return GNode::face_boundary(face_nodes);
    }
}

int& GCell::neib(int idx) {
    assert(idx < m_neibs.size());
    return m_neibs[idx];
}

int GCell::neib(int idx) const {
    assert(idx < m_neibs.size());
    return m_neibs[idx];
}


int Grid::n_cells() const {
    return m_cells.size();
}

void Grid::reserve_nodes(int size) {
    m_nodes.reserve(size);
}

void Grid::reserve_cells(int size) {
    m_cells.reserve(size);
}

void Grid::operator +=(GNode::Ref node) {
    m_nodes.emplace_back(node);
}

void Grid::operator +=(const GCell& cell) {
    m_cells.push_back(cell);
}

GNode::Ptr Grid::node(int idx) {
    assert(idx < m_nodes.size());
    return m_nodes[idx];
}

GNode::Ref Grid::node(int idx) const {
    assert(idx < m_nodes.size());
    return m_nodes[idx];
}

GCell& Grid::cell(int idx) {
    assert(idx < m_cells.size());
    return m_cells[idx];
}

const GCell& Grid::cell(int idx) const {
    assert(idx < m_cells.size());
    return m_cells[idx];
}

void Grid::setup_adjacency() {
    // Находим для вершин смежные ячейки
    for (auto& cell: m_cells) {
        cell.add_self_to_nodes();
    }

    // Находим соседство через грани
    for (auto& cell: m_cells) {
        cell.find_neibs();
    }
}

AmrCell Grid::amr_cell(int idx) const {
    auto gcell = cell(idx);

    assert(gcell.index == idx);

    if (gcell.type() == CellType::QUAD) {
        auto v1 = m_nodes[gcell.node(0).index];
        auto v2 = m_nodes[gcell.node(1).index];
        auto v3 = m_nodes[gcell.node(2).index];
        auto v4 = m_nodes[gcell.node(3).index];

        ShortList2D vlist = { v1->v, v2->v, v3->v, v4->v };

        AmrCell cell(vlist);

        cell.visualize();

        cell.rank  = mpi::rank();
        cell.index = idx;

        // Данные AMR
        cell.b_idx = idx;
        cell.z_idx = 0;
        cell.next  = 0;
        cell.level = 0;
        cell.flag  = 0;

        cell.faces[Side::L].boundary = gcell.boundary({v1, v3});
        cell.faces[Side::R].boundary = gcell.boundary({v2, v4});
        cell.faces[Side::B].boundary = gcell.boundary({v1, v2});
        cell.faces[Side::T].boundary = gcell.boundary({v3, v4});

        cell.faces[Side::L].adjacent.rank = 0;
        cell.faces[Side::R].adjacent.rank = 0;
        cell.faces[Side::B].adjacent.rank = 0;
        cell.faces[Side::T].adjacent.rank = 0;

        cell.faces[Side::L].adjacent.ghost = -1;
        cell.faces[Side::R].adjacent.ghost = -1;
        cell.faces[Side::B].adjacent.ghost = -1;
        cell.faces[Side::T].adjacent.ghost = -1;

        cell.faces[Side::L].adjacent.index = gcell.adjacent({v1, v3});
        cell.faces[Side::R].adjacent.index = gcell.adjacent({v2, v4});
        cell.faces[Side::B].adjacent.index = gcell.adjacent({v1, v2});
        cell.faces[Side::T].adjacent.index = gcell.adjacent({v3, v4});

        return cell;
    }
    else {
        throw std::runtime_error("Can't create AmrCell");
    }
}

} // namespace zephyr::geom

