#include <cassert>

#include <zephyr/utils/mpi.h>
#include <zephyr/geom/grid.h>

#include <zephyr/mesh/primitives/side.h>
#include <zephyr/mesh/primitives/mov_node.h>
#include <zephyr/mesh/primitives/amr_cell.h>
#include <zephyr/mesh/primitives/mov_cell.h>

using zephyr::mesh::Side3D;
using zephyr::mesh::AmrCell;
using zephyr::mesh::MovNode;
using zephyr::mesh::MovCell;

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

void GNode::set_boundaries(std::set<Boundary> flags) {
    m_bounds = flags;
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
        for (auto& node: vs) {
            if (!node->m_bounds.empty()) {
                return *node->m_bounds.begin();
            }
        }
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

MovNode GNode::bnode() const {
    MovNode mov_node(0, index);
    if (m_bounds.empty()) {
        mov_node.boundary = *m_bounds.begin();
    }
    mov_node.coords = this->v;

    int counter = 0;
    for (int idx: m_neibs) {
        mov_node.neibs[counter].rank = 0;
        mov_node.neibs[counter].index = idx;
        mov_node.neibs[counter].alien = -1;
        ++counter;
    }
    return mov_node;
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
    Vector3d a = nodes[1]->v - nodes[0]->v;
    Vector3d b = nodes[2]->v - nodes[1]->v;
    Vector3d c = nodes[3]->v - nodes[2]->v;
    Vector3d d = nodes[0]->v - nodes[3]->v;

    if (a.cross(b).dot(c.cross(d)) > 0.0) {
        // Верный порядок по часовой
        cell.m_nodes = {nodes[0], nodes[1], nodes[2], nodes[3]};
    } else {
        // Z-порядок обхода (меняем)
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

GCell GCell::polygon(std::vector<GNode::Ptr> nodes) {
    GCell cell(CellType::POLYGON);

    int n_nodes = nodes.size();

    Vector3d c = Vector3d::Zero();
    for (auto& n: nodes) {
        c += n->v;
    }
    c /= n_nodes;

    std::sort(nodes.begin(), nodes.end(),
              [&c](GNode::Ref a, GNode::Ref b) -> bool {
                  double phi_a = std::atan2(a->v.y() - c.y(), a->v.x() - c.x());
                  double phi_b = std::atan2(b->v.y() - c.y(), b->v.x() - c.x());
                  return phi_a < phi_b;
              });

    cell.m_nodes = nodes;

    cell.m_faces.resize(n_nodes);
    for (int i = 0; i < n_nodes; ++i) {
        cell.m_faces[i] = {i, (i + 1) % n_nodes};
    }

    cell.m_neibs.resize(n_nodes, -1);
    return cell;
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

GCell GCell::hexagedron(const std::array<GNode::Ptr, 8>& nodes) {
    GCell cell(CellType::HEXAHEDRON);

    // Проверка нумерации вершин в m_nodes
    Vector3d a = nodes[1]->v - nodes[0]->v;
    Vector3d b = nodes[2]->v - nodes[1]->v;
    Vector3d c = nodes[3]->v - nodes[2]->v;
    Vector3d d = nodes[0]->v - nodes[3]->v;

    if (a.cross(b).dot(c.cross(d)) > 0.0) {
        // Правильный порядок
        cell.m_nodes = {nodes[0], nodes[1], nodes[2], nodes[3],
                        nodes[4], nodes[5], nodes[6], nodes[7],};
    } else {
        // Z-порядок обхода (меняем)
        cell.m_nodes = {nodes[0], nodes[1], nodes[3], nodes[2],
                        nodes[4], nodes[5], nodes[7], nodes[6]};
    }
    cell.m_faces = {
            {0, 4, 7, 3},
            {1, 5, 6, 2},
            {0, 1, 5, 4},
            {3, 2, 6, 7},
            {0, 1, 2, 3},
            {4, 5, 6, 7}
    };
    cell.m_neibs = {-1, -1, -1, -1, -1, -1};
    return cell;
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

Grid::Grid() : structured(false) {

}

void Grid::set_axial(bool axial) {
    m_axial = axial;
}

int Grid::n_nodes() const {
    return m_nodes.size();
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
        auto v3 = m_nodes[gcell.node(3).index];
        auto v4 = m_nodes[gcell.node(2).index];

        Quad vlist = {v1->v, v2->v, v3->v, v4->v};

        AmrCell cell(vlist, m_axial);

        cell.rank  = mpi::rank();
        cell.index = idx;

        // Данные AMR
        cell.b_idx = idx;
        cell.z_idx = 0;
        cell.next  = 0;
        cell.level = 0;
        cell.flag  = 0;

        cell.faces[Side3D::L].boundary = gcell.boundary({v1, v3});
        cell.faces[Side3D::R].boundary = gcell.boundary({v2, v4});
        cell.faces[Side3D::B].boundary = gcell.boundary({v1, v2});
        cell.faces[Side3D::T].boundary = gcell.boundary({v3, v4});

        cell.faces[Side3D::L].adjacent.rank = 0;
        cell.faces[Side3D::R].adjacent.rank = 0;
        cell.faces[Side3D::B].adjacent.rank = 0;
        cell.faces[Side3D::T].adjacent.rank = 0;

        cell.faces[Side3D::L].adjacent.alien = -1;
        cell.faces[Side3D::R].adjacent.alien = -1;
        cell.faces[Side3D::B].adjacent.alien = -1;
        cell.faces[Side3D::T].adjacent.alien = -1;

        cell.faces[Side3D::L].adjacent.index = gcell.adjacent({v1, v3});
        cell.faces[Side3D::R].adjacent.index = gcell.adjacent({v2, v4});
        cell.faces[Side3D::B].adjacent.index = gcell.adjacent({v1, v2});
        cell.faces[Side3D::T].adjacent.index = gcell.adjacent({v3, v4});

        return cell;
    }
    else if (gcell.type() == CellType::HEXAHEDRON) {
        auto v0 = m_nodes[gcell.node(0).index];
        auto v1 = m_nodes[gcell.node(1).index];
        auto v2 = m_nodes[gcell.node(3).index];
        auto v3 = m_nodes[gcell.node(2).index];
        auto v4 = m_nodes[gcell.node(4).index];
        auto v5 = m_nodes[gcell.node(5).index];
        auto v6 = m_nodes[gcell.node(7).index];
        auto v7 = m_nodes[gcell.node(6).index];

        Cube vlist = {v0->v, v1->v, v2->v, v3->v,
                      v4->v, v5->v, v6->v, v7->v };

        AmrCell cell(vlist);

        cell.rank  = mpi::rank();
        cell.index = idx;

        // Данные AMR
        cell.b_idx = idx;
        cell.z_idx = 0;
        cell.next  = 0;
        cell.level = 0;
        cell.flag  = 0;

        cell.faces[Side3D::L].boundary = gcell.boundary({v0, v4, v2, v6});
        cell.faces[Side3D::R].boundary = gcell.boundary({v5, v1, v3, v7});
        cell.faces[Side3D::B].boundary = gcell.boundary({v0, v1, v4, v5});
        cell.faces[Side3D::T].boundary = gcell.boundary({v6, v7, v2, v3});
        cell.faces[Side3D::X].boundary = gcell.boundary({v0, v1, v2, v3});
        cell.faces[Side3D::F].boundary = gcell.boundary({v4, v5, v6, v7});

        cell.faces[Side3D::L].adjacent.rank = 0;
        cell.faces[Side3D::R].adjacent.rank = 0;
        cell.faces[Side3D::B].adjacent.rank = 0;
        cell.faces[Side3D::T].adjacent.rank = 0;
        cell.faces[Side3D::X].adjacent.rank = 0;
        cell.faces[Side3D::F].adjacent.rank = 0;

        cell.faces[Side3D::L].adjacent.alien = -1;
        cell.faces[Side3D::R].adjacent.alien = -1;
        cell.faces[Side3D::B].adjacent.alien = -1;
        cell.faces[Side3D::T].adjacent.alien = -1;
        cell.faces[Side3D::X].adjacent.alien = -1;
        cell.faces[Side3D::F].adjacent.alien = -1;

        cell.faces[Side3D::L].adjacent.index = gcell.adjacent({v0, v4, v2, v6});
        cell.faces[Side3D::R].adjacent.index = gcell.adjacent({v5, v1, v3, v7});
        cell.faces[Side3D::B].adjacent.index = gcell.adjacent({v0, v1, v4, v5});
        cell.faces[Side3D::T].adjacent.index = gcell.adjacent({v6, v7, v2, v3});
        cell.faces[Side3D::X].adjacent.index = gcell.adjacent({v0, v1, v2, v3});
        cell.faces[Side3D::F].adjacent.index = gcell.adjacent({v4, v5, v6, v7});

        return cell;
    }
    else if (gcell.type() == CellType::POLYGON) {
        int n_nodes = gcell.n_nodes();

        std::vector<GNode::Ptr> nodes(n_nodes, nullptr);

        Polygon poly(n_nodes);
        for (int i = 0; i < n_nodes; ++i) {
            nodes[i] = m_nodes[gcell.node(i).index];
            poly.set(i, nodes[i]->v);
        }

        AmrCell cell(poly);

        cell.rank  = mpi::rank();
        cell.index = idx;

        // Данные AMR
        cell.b_idx = idx;
        cell.z_idx = 0;
        cell.next  = 0;
        cell.level = 0;
        cell.flag  = 0;

        for (int i = 0; i < n_nodes; ++i) {
            GNode::Ref v1 = nodes[cell.faces[i].vertices[0]];
            GNode::Ref v2 = nodes[cell.faces[i].vertices[1]];

            cell.faces[i].boundary = gcell.boundary({v1, v2});
            cell.faces[i].adjacent.rank = 0;
            cell.faces[i].adjacent.alien = -1;
            cell.faces[i].adjacent.index = gcell.adjacent({v1, v2});
        }

        return cell;
    }
    else {
        throw std::runtime_error("Can't create AmrCell");
    }
}

MovCell Grid::mov_cell(int idx) const {
    auto gcell = cell(idx);

    assert(gcell.index == idx);

    if (gcell.type() == CellType::QUAD ||
        gcell.type() == CellType::POLYGON) {

        int n_nodes = gcell.n_nodes();

        std::vector<GNode::Ptr> nodes(n_nodes, nullptr);

        geom::Polygon poly(n_nodes);
        for (int i = 0; i < n_nodes; ++i) {
            nodes[i] = m_nodes[gcell.node(i).index];
            poly.set(i, nodes[i]->v);
        }

        MovCell cell(poly);
        for (int i = 0; i < n_nodes; ++i) {
            cell.nodes[i] = nodes[i]->index;
        }

        cell.rank  = mpi::rank();
        cell.index = idx;
        cell.next  = 0;

        for (int i = 0; i < n_nodes; ++i) {
            GNode::Ref v1 = nodes[cell.faces[i].vertices[0]];
            GNode::Ref v2 = nodes[cell.faces[i].vertices[1]];

            cell.faces[i].boundary = gcell.boundary({v1, v2});
            cell.faces[i].adjacent.rank = 0;
            cell.faces[i].adjacent.alien = -1;
            cell.faces[i].adjacent.index = gcell.adjacent({v1, v2});
        }

        return cell;
    }
    else {
        throw std::runtime_error("Can't create MovCell");
    }
}

void Grid::assume_structured(int nx, int ny, int nz) {
    if (nx <= 0 && ny <= 0 && nz <= 0) {
        throw std::runtime_error("Grid::assume structured nx | ny | nz <= 0");
    }

    if (m_cells.size() != nx * ny * nz) {
        throw std::runtime_error("Grid::assume structured nx * ny * nz != n_cells");
    }

    structured = true;
    m_nx = nx;
    m_ny = ny;
    m_nz = nz;
}

bool Grid::is_structured() const {
    return structured;
}

int Grid::nx() const {
    return structured ? m_nx : int(m_cells.size());
}

int Grid::ny() const {
    return structured ? m_ny : 1;
}

int Grid::nz() const {
    return structured ? m_nz : 1;
}

} // namespace zephyr::geom

