#include <limits>
#include <stdexcept>
#include <utility>
#include <format>
#include <ranges>
#include <unordered_set>

#include <zephyr/utils/threads.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/side.h>
#include <zephyr/geom/indexing.h>

namespace zephyr::geom {

using utils::threads;

void Grid::move(const Vector3d& shift) {
    require_editable_();
    threads::for_each(
        m_nodes.begin(), m_nodes.end(),
        [shift](Node& node) {
           node.pos += shift;
        });
}

void Grid::scale(double q) {
    require_editable_();
    threads::for_each(
        m_nodes.begin(), m_nodes.end(),
        [q](Node& node) {
           node.pos *= q;
        });
}

void Grid::rotate(double phi) {
    require_editable_();
    double sin = std::sin(phi);
    double cos = std::cos(phi);
    threads::for_each(
        m_nodes.begin(), m_nodes.end(),
        [sin, cos](Node& node) {
            double x = cos * node.pos.x() - sin * node.pos.y();
            double y = sin * node.pos.x() + cos * node.pos.y();
            node.pos.x() = x;
            node.pos.y() = y;
        });
}

void Grid::rotate(const Matrix3d& R) {
    require_editable_();
    threads::for_each(
        m_nodes.begin(), m_nodes.end(),
        [R](Node& node) {
           node.pos.applyOnTheLeft(R);
        });
}

void Grid::transform(std::function<Vector3d(const Vector3d&)>& func) {
    require_editable_();
    threads::for_each(
        m_nodes.begin(), m_nodes.end(),
        [&func](Node& node) {
           node.pos = func(node.pos);
        });
}

std::vector<Node::Ptr> Grid::create_central_nodes_() {
    std::vector<Node::Ptr> central(m_cells.size());
    for (id_t i = 0; i < m_cells.size(); ++i) {
        Vector3d c = m_cells[i].center(m_nodes);
        central[i] = Node::create(c);
    }
    for (const auto& node: central) {
        add_node(node);
    }
    return central;
}

std::unordered_map<FaceKey, Node::Ptr> Grid::create_face_nodes_() {
    std::unordered_map<FaceKey, Node::Ptr> face_nodes;
    for (auto & cell: m_cells) {
        if (!cell.has_faces()) {
            cell.init_faces();
        }
        for (int iface = 0; iface < cell.n_faces(); ++iface) {
            const Face& face = cell.get_face(iface);
            Vector3d fc = cell.face_center(m_nodes, iface);
            FaceKey face_key{cell, face};
            face_nodes[face_key] = Node::create(fc);
        }
    }
    for (const auto& val: face_nodes | std::views::values) {
        add_node(val);
    }
    return face_nodes;
}

std::unordered_map<EdgeKey, Node::Ptr> Grid::create_edge_nodes_() {
    if (m_dim != 3 && m_type != Type::QUAD) {
        throw std::runtime_error("Grid::create_edge_nodes_: Only for HEX cells");
    }

    std::unordered_map<EdgeKey, Node::Ptr> edge_nodes;
    edge_nodes.reserve(4 * m_cells.size()); // Для сетки из HEX ячеек
    for (auto& cell: m_cells) {
        if (!cell.has_faces()) {
            cell.init_faces();
        }
        for (const auto& face: cell.faces()) {
            for (int i = 0; i < face.n_nodes(); ++i) {
                int j = (i + 1) % face.n_nodes();
                id_t nid1 = cell.node_idx(face.node_idx(i));
                id_t nid2 = cell.node_idx(face.node_idx(j));
                EdgeKey edge{nid1, nid2};
                if (!edge_nodes.contains(edge)) {
                    Vector3d v1 = m_nodes[nid1].pos;
                    Vector3d v2 = m_nodes[nid2].pos;
                    edge_nodes[edge] = Node::create(0.5 * (v1 + v2));
                }
            }
        }
    }
    for (const auto& val: edge_nodes | std::views::values) {
        add_node(val);
    }
    return edge_nodes;
}

void Grid::mirror_(int axis) {
    require_editable_();
    if (m_nodes.empty() || m_cells.empty()) {
        std::cerr << "Grid::mirror: empty grid\n";
        return;
    }

    constexpr double eps = 1.0e-9;

    // Определить полуплоскость, в которой лежат узлы
    double sign = NAN;
    for (const auto& node: m_nodes) {
        if (std::abs(node.pos[axis]) > eps) {
            sign = node.pos[axis] < 0.0 ? -1.0 : +1.0;
            break;
        }
    }

    // Узлы на границе
    std::unordered_set<id_t> border_nodes;
    for (const auto& node: m_nodes) {
        if (std::abs(node.pos[axis]) < eps) {
            border_nodes.insert(node.id());
        }
        if (sign * node.pos[axis] < -eps) {
            throw std::runtime_error(std::format("Grid::mirror: node[axis] = {:.3e} is out of bound", node.pos[axis]));
        }
    }

    // Добавить отраженные вершины
    const id_t n_nodes_prev = m_nodes.size();
    std::vector<Node::Ptr> mirror_nodes(m_nodes.size());
    for (id_t i = 0; i < m_nodes.size(); ++i) {
        if (border_nodes.contains(i)) {
            mirror_nodes[i] = Node::create(m_nodes[i].pos);
            mirror_nodes[i]->m_id = m_nodes[i].id();
        }
        else {
            mirror_nodes[i] = Node::create(m_nodes[i].pos);
            mirror_nodes[i]->pos[axis] *= -1.0;
        }
    }

    m_nodes.reserve(2 * n_nodes_prev);
    for (const auto& node: mirror_nodes) {
        add_node(node);
    }

    // Создать копии существующих ячеек с индексами новых вершин
    const id_t n_cells_prev = m_cells.size();
    m_cells.reserve(2 * n_cells_prev);
    for (const auto& cell: m_cells) {
        m_cells.emplace_back(cell);
    }
    for (id_t ic = n_cells_prev; ic < m_cells.size(); ++ic) {
        std::vector<id_t> new_nodes = m_cells[ic].nodes();
        for (id_t& j: new_nodes) {
            j = mirror_nodes[j]->id();
        }
        m_cells[ic].replace_nodes(std::move(new_nodes));
        m_cells[ic].mirror();
    }
}

void Grid::mirror_x() {
    mirror_(0);
}

void Grid::mirror_y() {
    mirror_(1);
}

void Grid::mirror_z() {
    mirror_(2);
}

void Grid::triangulation(int mode) {
    require_editable_();
    if (m_cells.empty() || m_nodes.empty()) {
        std::cerr << "Grid::triangulation: empty grid\n";
        return;
    }
    if (m_type == Type::TRI) {
        std::cerr << "Grid::triangulation: grid is already triangulated\n";
        return;
    }
    if (m_type != Type::QUAD) {
        std::cerr << "Grid::triangulation: can only triangulate QUAD grid\n";
        return;
    }
    if (mode != 1 && mode != 2) {
        std::cerr << "Grid::triangulation: unknown triangulation mode " << mode << "\n";
    }

    if (dimension() == 2) {
        if (mode == 1) {
            triangulation_quad_delaunay_();
        }
        else {
            triangulation_quad_symmetry_();
        }
    }
    else {
        triangulation_hex_symmetry_();
    }
}

// Функция для вычисления определителя 4x4 для проверки, лежит ли d внутри окружности abc
// Определяет, лежит ли точка d внутри окружности, проходящей через a, b, c
double circle_determinant(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& d) {
    double a2 = a.squaredNorm();
    double b2 = b.squaredNorm();
    double c2 = c.squaredNorm();
    double d2 = d.squaredNorm();

    return a.x() * (b.y() * (c2 - d2) - c.y() * (b2 - d2) + d.y() * (b2 - c2)) -
           a.y() * (b.x() * (c2 - d2) - c.x() * (b2 - d2) + d.x() * (b2 - c2)) +
           a2 * (b.x() * c.y() - c.x() * b.y() - b.x() * d.y() + d.x() * b.y() + c.x() * d.y() - d.x() * c.y()) +
           b2 * (c.x() * d.y() - d.x() * c.y()) +
           c2 * (d.x() * a.y() - a.x() * d.y()) +
           d2 * (a.x() * b.y() - b.x() * a.y());
}

// Проверяет, является ли диагональ v1v3 триангуляцией Делоне
bool is_delaunay_diagonal_v1v3(const std::array<Vector3d, 4>& vs) {
    return circle_determinant(vs[0], vs[1], vs[2], vs[3]) < 0.0;
}

// Выбрать диагональ, которая достраивает триангуляцию Делоне.
// Возвращает 0 — диагональ v1v3, 1 — диагональ v2v4
int choose_delaunay_diagonal(const std::array<Vector3d, 4>& vs) {
    if (is_delaunay_diagonal_v1v3(vs)) {
        return 0; // Диагональ v1v3
    } else {
        return 1; // Диагональ v2v4
    }
}

void Grid::triangulation_quad_delaunay_() {
    const std::vector<Cell> quads = std::move(m_cells);

    id_t n_cells = quads.size();
    m_cells.clear();
    m_cells.reserve(2 * n_cells);
    for (id_t ic = 0; ic < n_cells; ++ic) {
        const auto& node_ids = quads[ic].nodes();
        auto face_bc = quads[ic].faces_bc<4>();
        std::array<Vector3d, 4> node_pos;
        for (int k = 0; k < 4; ++k) {
            node_pos[k] = m_nodes[node_ids[k]].pos;
        }

        using indexing::quad::vs;

        if (choose_delaunay_diagonal(node_pos) == 0) {
            Cell c1(CellType::TRIANGLE, {node_ids[vs<0, 0>()], node_ids[vs<1, 0>()], node_ids[vs<1, 1>()]});
            c1.set_face_bc({face_bc[Side2D::B], face_bc[Side2D::R], Boundary::INNER});

            Cell c2(CellType::TRIANGLE, {node_ids[vs<0, 0>()], node_ids[vs<1, 1>()], node_ids[vs<0, 1>()]});
            c2.set_face_bc({Boundary::INNER, face_bc[Side2D::T], face_bc[Side2D::L]});

            m_cells.emplace_back(std::move(c1));
            m_cells.emplace_back(std::move(c2));
        }
        else {
            Cell c1(CellType::TRIANGLE, {node_ids[vs<0, 0>()], node_ids[vs<1, 0>()], node_ids[vs<0, 1>()]});
            c1.set_face_bc({face_bc[Side2D::B], Boundary::INNER, face_bc[Side2D::L]});

            Cell c2(CellType::TRIANGLE, {node_ids[vs<1, 0>()], node_ids[vs<1, 1>()], node_ids[vs<0, 1>()]});
            c2.set_face_bc({face_bc[Side2D::R], face_bc[Side2D::T], Boundary::INNER});

            m_cells.emplace_back(std::move(c1));
            m_cells.emplace_back(std::move(c2));
        }
    }
    m_type = Type::TRI;
}

void Grid::triangulation_quad_symmetry_() {
    // Создать центральные вершины
    const auto central = create_central_nodes_();

    // Переместить исходные ячейки
    const std::vector<Cell> quads = std::move(m_cells);

    id_t n_cells = quads.size();
    m_cells.clear();
    m_cells.reserve(4 * n_cells);
    for (id_t ic = 0; ic < n_cells; ++ic) {
        id_t central_id = central[ic]->id();

        const auto& node_ids = quads[ic].nodes();
        auto face_bc = quads[ic].faces_bc<4>();

        using indexing::quad::vs;
        Cell c1(CellType::TRIANGLE, {node_ids[vs<0, 1>()], node_ids[vs<0, 0>()], central_id});
        c1.set_face_bc({face_bc[Side2D::L], Boundary::INNER, Boundary::INNER});

        Cell c2(CellType::TRIANGLE, {node_ids[vs<1, 0>()], node_ids[vs<1, 1>()], central_id});
        c2.set_face_bc({face_bc[Side2D::R], Boundary::INNER, Boundary::INNER});

        Cell c3(CellType::TRIANGLE, {node_ids[vs<0, 0>()], node_ids[vs<1, 0>()], central_id});
        c3.set_face_bc({face_bc[Side2D::B], Boundary::INNER, Boundary::INNER});

        Cell c4(CellType::TRIANGLE, {node_ids[vs<1, 1>()], node_ids[vs<0, 1>()], central_id});
        c4.set_face_bc({face_bc[Side2D::T], Boundary::INNER, Boundary::INNER});

        m_cells.emplace_back(std::move(c1));
        m_cells.emplace_back(std::move(c2));
        m_cells.emplace_back(std::move(c3));
        m_cells.emplace_back(std::move(c4));
    }
    m_type = Type::TRI;
}

void Grid::triangulation_hex_symmetry_() {
    require_editable_();
    if (m_cells.empty() || m_nodes.empty()) {
        std::cerr << "Grid::triangulation: empty grid\n";
        return;
    }
    if (m_type != Type::QUAD) {
        std::cerr << "Grid::triangulation: can only pyramidize QUAD grid\n";
        return;
    }
    if (m_dim != 3) {
        throw std::runtime_error("Unreachable");
    }

    // Создать центральные вершины
    const auto central = create_central_nodes_();

    // Создать вершины граней
    auto face_unique_nodes = create_face_nodes_();

    // Переместить исходные ячейки
    const std::vector<Cell> cubes = std::move(m_cells);

    id_t n_cells = cubes.size();
    m_cells.clear();
    m_cells.reserve(6 * n_cells);
    for (id_t i = 0; i < n_cells; ++i) {
        id_t central_id = central[i]->id();

        const auto& node_ids = cubes[i].nodes();
        auto face_bc = cubes[i].faces_bc<6>();

        std::array<id_t, 6> face_ids;
        for (int k = 0; k < 6; ++k) {
            FaceKey key{cubes[i], cubes[i].get_face(k)};
            face_ids[k] = face_unique_nodes[key]->id();
        }

        for (int side = 0; side < 6; ++side) {
            // У hex грани обходятся нормалью наружу
            using indexing::hex::face_nodes;

            Cell c1(CellType::TETRA, {face_ids[side], node_ids[face_nodes[side][1]], node_ids[face_nodes[side][0]], central_id});
            c1.set_face_bc({face_bc[side], Boundary::INNER, Boundary::INNER, Boundary::INNER});

            Cell c2(CellType::TETRA, {face_ids[side], node_ids[face_nodes[side][2]], node_ids[face_nodes[side][1]], central_id});
            c2.set_face_bc({face_bc[side], Boundary::INNER, Boundary::INNER, Boundary::INNER});

            Cell c3(CellType::TETRA, {face_ids[side], node_ids[face_nodes[side][3]], node_ids[face_nodes[side][2]], central_id});
            c3.set_face_bc({face_bc[side], Boundary::INNER, Boundary::INNER, Boundary::INNER});

            Cell c4(CellType::TETRA, {face_ids[side], node_ids[face_nodes[side][0]], node_ids[face_nodes[side][3]], central_id});
            c4.set_face_bc({face_bc[side], Boundary::INNER, Boundary::INNER, Boundary::INNER});

            m_cells.emplace_back(std::move(c1));
            m_cells.emplace_back(std::move(c2));
            m_cells.emplace_back(std::move(c3));
            m_cells.emplace_back(std::move(c4));
        }
    }
    m_type = Type::TRI;
}

void Grid::pyramidize() {
    require_editable_();
    if (m_cells.empty() || m_nodes.empty()) {
        std::cerr << "Grid::pyramidize: empty grid\n";
        return;
    }
    if (m_type != Type::QUAD) {
        std::cerr << "Grid::pyramidize: can only pyramidize QUAD grid\n";
        return;
    }
    if (m_dim == 2) {
        triangulation(2);
        return;
    }

    // Создать центральные вершины
    const auto central = create_central_nodes_();

    // Переместить исходные ячейки
    const std::vector<Cell> cubes = std::move(m_cells);

    id_t n_cells = cubes.size();
    m_cells.clear();
    m_cells.reserve(6 * n_cells);
    for (id_t i = 0; i < n_cells; ++i) {
        id_t central_id = central[i]->id();

        const auto& base_ids = cubes[i].nodes();
        auto face_bc = cubes[i].faces_bc<6>();

        for (int side = 0; side < 6; ++side) {
            // У hex грани обходятся нормалью наружу, у пирамиды
            // узлы основания перечисляются нормалью внутрь
            using indexing::hex::face_nodes;

            Cell cell(CellType::PYRAMID, {
                base_ids[face_nodes[side][3]],
                base_ids[face_nodes[side][2]],
                base_ids[face_nodes[side][1]],
                base_ids[face_nodes[side][0]],
                central_id});

            cell.set_face_bc({face_bc[side],
                Boundary::INNER, Boundary::INNER,
                Boundary::INNER, Boundary::INNER
            });

            m_cells.emplace_back(std::move(cell));
        }
    }
    m_type = Type::POLY;
}

void Grid::extrude(const Vector3d& p, int N, Boundary side1, Boundary side2) {
    require_editable_();
    if (m_cells.empty() || m_nodes.empty()) {
        std::cerr << "Grid::extrude: empty grid\n";
        return;
    }
    if (m_dim != 2) {
        std::cerr << "Grid::extrude: can't extrude from 3D grid\n";
        return;
    }
    if (m_type == Type::AMR) {
        std::cerr << "Grid::extrude: can't extrude AMR grid\n";
        return;
    }
    if (p.z() < 0.0) {
        std::cerr << "Grid::extrude: negative extrude not implemented\n";
        return;
    }

    m_dim = 3;

    // Создать новые вершины
    id_t n_nodes = m_nodes.size();
    m_nodes.reserve((N + 1) * n_nodes);
    for (id_t k = 1; k <= N; ++k) {
        for (id_t j = 0; j < n_nodes; ++j) {
            auto node = Node::create(m_nodes[j].pos + (k * p) / N);
            add_node(node);
        }
    }

    const std::vector<Cell> grid = std::move(m_cells);

    id_t n_cells = grid.size();
    m_cells.reserve(N * n_cells);
    for (id_t ic = 0; ic < n_cells; ++ic) {
        CellType type = grid[ic].type();
        const auto& base_ids = grid[ic].nodes();
        std::vector base_bc(grid[ic].n_faces(), Boundary::UNDEFINED);
        for (int k = 0; k < base_bc.size(); ++k) {
            base_bc[k] = grid[ic].get_face(k).bc();
        }

        for (id_t k = 0; k < N; ++k) {
            // MeshType - TRI, QUAD, POLY, possible CellType is the same
            if (type == CellType::TRIANGLE) {
                // Преобразуем в WEDGE (треугольная призма)
                using namespace indexing;
                Cell cell(CellType::WEDGE, {
                    base_ids[0] + k * n_nodes,
                    base_ids[2] + k * n_nodes,
                    base_ids[1] + k * n_nodes,
                    base_ids[0] + (k + 1) * n_nodes,
                    base_ids[2] + (k + 1) * n_nodes,
                    base_ids[1] + (k + 1) * n_nodes
                });
                cell.set_face_bc({
                    base_bc[2],
                    base_bc[1],
                    base_bc[0],
                    k == 0 ? side1 : Boundary::INNER,
                    (k == N - 1) ? side2 : Boundary::INNER
                });
                m_cells.emplace_back(std::move(cell));
                m_type = Type::POLY;
            }
            else if (type == CellType::QUAD) {
                // Преобразуем в HEX
                using namespace indexing;
                Cell cell(CellType::HEXAHEDRON, {
                    base_ids[quad::vs<0,0>()] + k * n_nodes,
                    base_ids[quad::vs<1,0>()] + k * n_nodes,
                    base_ids[quad::vs<1,1>()] + k * n_nodes,
                    base_ids[quad::vs<0,1>()] + k * n_nodes,
                    base_ids[quad::vs<0,0>()] + (k + 1) * n_nodes,
                    base_ids[quad::vs<1,0>()] + (k + 1) * n_nodes,
                    base_ids[quad::vs<1,1>()] + (k + 1) * n_nodes,
                    base_ids[quad::vs<0,1>()] + (k + 1) * n_nodes
                });
                cell.set_face_bc({
                    base_bc[Side2D::L],
                    base_bc[Side2D::R],
                    base_bc[Side2D::B],
                    base_bc[Side2D::T],
                    k == 0 ? side1 : Boundary::INNER,
                    (k == N - 1) ? side2 : Boundary::INNER
                });
                m_cells.emplace_back(std::move(cell));
                if (m_type == Type::TRI) { m_type = Type::POLY; }
            }
            else if (type == CellType::POLYGON) {
                z_assert(base_bc.size() == base_ids.size(), "Wrong polygon");
                int poly_size = base_ids.size();

                // Вершины многогранника
                std::vector<id_t> node_ids(2 * poly_size);
                for (int i = 0; i < poly_size; ++i) {
                    node_ids[i] = base_ids[i] + k * n_nodes;
                    node_ids[i + poly_size] = base_ids[i] + (k + 1) * n_nodes;
                }

                Cell cell(CellType::POLYHEDRON, std::move(node_ids));

                // Массив граничных условий
                std::vector face_bc(poly_size + 2, Boundary::UNDEFINED);
                for (int i = 0; i < poly_size; ++i) {
                    face_bc[i] = base_bc[i];
                }
                face_bc[poly_size + 0] = k == 0 ? side1 : Boundary::INNER;
                face_bc[poly_size + 1] = (k == N - 1) ? side2 : Boundary::INNER;

                std::vector<std::vector<int>> face_indices(poly_size + 2);
                for (int i = 0; i < poly_size; ++i) {
                    int j = (i + 1) % poly_size;
                    face_indices[i].reserve(4);
                    face_indices[i].emplace_back(i);
                    face_indices[i].emplace_back(j);
                    face_indices[i].emplace_back(j + poly_size);
                    face_indices[i].emplace_back(i + poly_size);
                }
                face_indices[poly_size + 0].reserve(poly_size);
                for (int i = poly_size - 1; i >= 0; --i) {
                    face_indices[poly_size + 0].emplace_back(i);
                }
                face_indices[poly_size + 1].reserve(poly_size);
                for (int i = 0; i < poly_size; ++i) {
                    face_indices[poly_size + 1].emplace_back(poly_size + i);
                }
                cell.set_faces(face_indices);
                cell.set_face_bc(face_bc);
                m_cells.emplace_back(std::move(cell));
                m_type = Type::POLYHEDRON;
            }
            else {
                throw std::runtime_error("Grid::extrude: unknown cell type");
            }
        }
    }
}

void Grid::make_amr_2D_() {
    // Создать центральные вершины
    const auto central = create_central_nodes_();

    // Создать вершины граней
    auto face_nodes = create_face_nodes_();

    // Сохраняем исходные ячейки
    const std::vector<Cell> quads = std::move(m_cells);

    // Преобразовать ячейки в AMR2D (заменяем полигоны на таблицы 3x3)
    m_cells.clear();
    for (id_t i = 0; i < quads.size(); ++i) {
        const Cell& quad = quads[i];
        FaceKey L{quad, quad.get_face(Side2D::L)};
        FaceKey R{quad, quad.get_face(Side2D::R)};
        FaceKey B{quad, quad.get_face(Side2D::B)};
        FaceKey T{quad, quad.get_face(Side2D::T)};

        z_assert(face_nodes.contains(L), "Cell face not found");
        z_assert(face_nodes.contains(R), "Cell face not found");
        z_assert(face_nodes.contains(B), "Cell face not found");
        z_assert(face_nodes.contains(T), "Cell face not found");

        id_t vLB = quad.node_idx(0);
        id_t vRB = quad.node_idx(1);
        id_t vLT = quad.node_idx(3);
        id_t vRT = quad.node_idx(2);
        id_t v_C = central[i]->id();
        id_t v_L = face_nodes[L]->id();
        id_t v_R = face_nodes[R]->id();
        id_t v_B = face_nodes[B]->id();
        id_t v_T = face_nodes[T]->id();

        Cell cell(CellType::AMR2D,{
            vLB, v_B, vRB,
            v_L, v_C, v_R,
            vLT, v_T, vRT
        });
        cell.init_faces();
        for (int k = 0; k < 4; ++k) {
            cell.set_bc(k, quad.get_face(k).bc());
        }
        m_cells.emplace_back(std::move(cell));
    }
    m_type = Type::AMR;
}

void Grid::make_amr_3D_() {
    // Создать центральные вершины
    const auto central = create_central_nodes_();

    // Создать вершины граней
    auto face_nodes = create_face_nodes_();

    // Создать узлы на рёбрах
    auto edge_nodes = create_edge_nodes_();

    // Сохраняем исходные ячейки
    const std::vector<Cell> quads = std::move(m_cells);

    // Преобразовать ячейки в AMR2D (заменяем полигоны на таблицы 3x3)
    m_cells.clear();
    for (id_t i = 0; i < quads.size(); ++i) {
        const Cell& quad = quads[i];
        FaceKey L{quad, quad.get_face(Side3D::L)};
        FaceKey R{quad, quad.get_face(Side3D::R)};
        FaceKey B{quad, quad.get_face(Side3D::B)};
        FaceKey T{quad, quad.get_face(Side3D::T)};
        FaceKey Z{quad, quad.get_face(Side3D::Z)};
        FaceKey F{quad, quad.get_face(Side3D::F)};

        z_assert(face_nodes.contains(L), "Cell face not found");
        z_assert(face_nodes.contains(R), "Cell face not found");
        z_assert(face_nodes.contains(B), "Cell face not found");
        z_assert(face_nodes.contains(T), "Cell face not found");
        z_assert(face_nodes.contains(Z), "Cell face not found");
        z_assert(face_nodes.contains(F), "Cell face not found");

        using namespace indexing;
        
        // Ребра вдоль Ox
        EdgeKey e_bz{quad.node_idx(hex::vs<0,0,0>()), quad.node_idx(hex::vs<1,0,0>())};
        EdgeKey e_bf{quad.node_idx(hex::vs<0,0,1>()), quad.node_idx(hex::vs<1,0,1>())};
        EdgeKey e_tz{quad.node_idx(hex::vs<0,1,0>()), quad.node_idx(hex::vs<1,1,0>())};
        EdgeKey e_tf{quad.node_idx(hex::vs<0,1,1>()), quad.node_idx(hex::vs<1,1,1>())};

        // Ребра вдоль Oy
        EdgeKey e_zl{quad.node_idx(hex::vs<0,0,0>()), quad.node_idx(hex::vs<0,1,0>())};
        EdgeKey e_zr{quad.node_idx(hex::vs<1,0,0>()), quad.node_idx(hex::vs<1,1,0>())};
        EdgeKey e_fl{quad.node_idx(hex::vs<0,0,1>()), quad.node_idx(hex::vs<0,1,1>())};
        EdgeKey e_fr{quad.node_idx(hex::vs<1,0,1>()), quad.node_idx(hex::vs<1,1,1>())};

        // Ребра вдоль Oz
        EdgeKey e_lb{quad.node_idx(hex::vs<0,0,0>()), quad.node_idx(hex::vs<0,0,1>())};
        EdgeKey e_lt{quad.node_idx(hex::vs<0,1,0>()), quad.node_idx(hex::vs<0,1,1>())};
        EdgeKey e_rb{quad.node_idx(hex::vs<1,0,0>()), quad.node_idx(hex::vs<1,0,1>())};
        EdgeKey e_rt{quad.node_idx(hex::vs<1,1,0>()), quad.node_idx(hex::vs<1,1,1>())};

        std::vector<id_t> nodes(27);
        // Угловые вершины
        nodes[amr::vs<-1, -1, -1>()] = quad.node_idx(hex::vs<0, 0, 0>());
        nodes[amr::vs<+1, -1, -1>()] = quad.node_idx(hex::vs<1, 0, 0>());
        nodes[amr::vs<-1, +1, -1>()] = quad.node_idx(hex::vs<0, 1, 0>());
        nodes[amr::vs<+1, +1, -1>()] = quad.node_idx(hex::vs<1, 1, 0>());
        nodes[amr::vs<-1, -1, +1>()] = quad.node_idx(hex::vs<0, 0, 1>());
        nodes[amr::vs<+1, -1, +1>()] = quad.node_idx(hex::vs<1, 0, 1>());
        nodes[amr::vs<-1, +1, +1>()] = quad.node_idx(hex::vs<0, 1, 1>());
        nodes[amr::vs<+1, +1, +1>()] = quad.node_idx(hex::vs<1, 1, 1>());

        // Центральная
        nodes[amr::vs<0, 0, 0>()] = central[i]->id();

        // Вершины на гранях
        nodes[amr::vs<-1, 0, 0>()] = face_nodes[L]->id();
        nodes[amr::vs<+1, 0, 0>()] = face_nodes[R]->id();
        nodes[amr::vs<0, -1, 0>()] = face_nodes[B]->id();
        nodes[amr::vs<0, +1, 0>()] = face_nodes[T]->id();
        nodes[amr::vs<0, 0, -1>()] = face_nodes[Z]->id();
        nodes[amr::vs<0, 0, +1>()] = face_nodes[F]->id();

        // Вершины на рёбрах
        nodes[amr::vs<0, -1, -1>()] = edge_nodes[e_bz]->id();
        nodes[amr::vs<0, -1, +1>()] = edge_nodes[e_bf]->id();
        nodes[amr::vs<0, +1, -1>()] = edge_nodes[e_tz]->id();
        nodes[amr::vs<0, +1, +1>()] = edge_nodes[e_tf]->id();

        // Ребра вдоль Oy
        nodes[amr::vs<-1, 0, -1>()] = edge_nodes[e_zl]->id();
        nodes[amr::vs<+1, 0, -1>()] = edge_nodes[e_zr]->id();
        nodes[amr::vs<-1, 0, +1>()] = edge_nodes[e_fl]->id();
        nodes[amr::vs<+1, 0, +1>()] = edge_nodes[e_fr]->id();

        // Ребра вдоль Oz
        nodes[amr::vs<-1, -1, 0>()] = edge_nodes[e_lb]->id();
        nodes[amr::vs<-1, +1, 0>()] = edge_nodes[e_lt]->id();
        nodes[amr::vs<+1, -1, 0>()] = edge_nodes[e_rb]->id();
        nodes[amr::vs<+1, +1, 0>()] = edge_nodes[e_rt]->id();

        Cell cell(CellType::AMR3D, std::move(nodes));
        cell.init_faces();
        for (int k = 0; k < 6; ++k) {
            cell.set_bc(k, quad.get_face(k).bc());
        }
        m_cells.emplace_back(std::move(cell));
    }
    m_type = Type::AMR;
}

void Grid::Grid::make_amr() {
    require_editable_();
    if (m_cells.empty() || m_nodes.empty()) {
        std::cerr << "Grid::make_amr: empty grid\n";
        return;
    }
    if (m_type != Type::QUAD) {
        throw std::runtime_error("Grid::make_amr: grid is not QUAD type");
    }
    if (dimension() == 2) {
        make_amr_2D_();
    }
    else {
        make_amr_3D_();
    }
}

} // namespace zephyr::geom