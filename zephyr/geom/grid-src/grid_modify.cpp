#include <limits>
#include <stdexcept>
#include <utility>
#include <format>
#include <ranges>
#include <unordered_set>

#include <zephyr/utils/threads.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/side.h>

namespace zephyr::geom {

using utils::threads;

void Grid::move(const Vector3d& shift) {
    require_editable_();
    if (!m_draft) {
        throw std::runtime_error("Grid::move: no draft data");
    }

    DraftData& draft = *m_draft;
    threads::for_each(
        draft.nodes.begin(), draft.nodes.end(),
        [shift](NodeInput::Ref node) {
           node->pos += shift;
        });
}

void Grid::scale(double q) {
    require_editable_();
    if (!m_draft) {
        throw std::runtime_error("Grid::scale: no draft data");
    }

    DraftData& draft = *m_draft;
    threads::for_each(
        draft.nodes.begin(), draft.nodes.end(),
        [q](NodeInput::Ref node) {
           node->pos *= q;
        });
}

void Grid::rotate(double phi) {
    require_editable_();
    if (!m_draft) {
        throw std::runtime_error("Grid::rotate: no draft data");
    }
    DraftData& draft = *m_draft;
    double sin = std::sin(phi);
    double cos = std::cos(phi);
    threads::for_each(
        draft.nodes.begin(), draft.nodes.end(),
        [sin, cos](NodeInput::Ref node) {
            double x = cos * node->pos.x() - sin * node->pos.y();
            double y = sin * node->pos.x() + cos * node->pos.y();
            node->pos.x() = x;
            node->pos.y() = y;
        });
}

void Grid::rotate(const Matrix3d& R) {
    require_editable_();
    if (!m_draft) {
        throw std::runtime_error("Grid::rotate: no draft data");
    }
    DraftData& draft = *m_draft;
    threads::for_each(
        draft.nodes.begin(), draft.nodes.end(),
        [R](NodeInput::Ref node) {
           node->pos.applyOnTheLeft(R);
        });
}

void Grid::transform(std::function<Vector3d(const Vector3d&)>& func) {
    require_editable_();
    if (!m_draft) {
        throw std::runtime_error("Grid::transform: no draft data");
    }
    DraftData& draft = *m_draft;
    threads::for_each(
        draft.nodes.begin(), draft.nodes.end(),
        [&func](NodeInput::Ref node) {
           node->pos = func(node->pos);
        });
}

void Grid::mirror_(int axis) {
    require_editable_();
    if (!m_draft) {
        throw std::runtime_error("Grid::mirror: no draft data");
    }
    DraftData& draft = *m_draft;

    constexpr double eps = 1.0e-9;

    if (draft.nodes.empty()) {
        std::cerr << "Grid::mirror: empty nodes, skip\n";
        return;
    }

    // Определить полуплоскость, в которой лежат узлы
    double sign = NAN;
    for (const auto node: draft.nodes) {
        if (std::abs(node->pos[axis]) > eps) {
            sign = node->pos[axis] < 0.0 ? -1.0 : +1.0;
            break;
        }
    }

    // Узлы на границе
    std::unordered_set<id_t> border_nodes;
    for (const auto node: draft.nodes) {
        if (std::abs(node->pos[axis]) < eps) {
            border_nodes.insert(node->m_builder.id);
        }
        if (sign * node->pos[axis] < -eps) {
            throw std::runtime_error(std::format("Grid::mirror: node[axis] = {:.3e} is out of bound", node->pos[axis]));
        }
    }

    // Добавить отраженные вершины
    const id_t n_nodes_prev = draft.nodes.size();
    std::vector<NodeInput::Ptr> mirror_nodes(draft.nodes.size());
    for (id_t i = 0; i < draft.nodes.size(); ++i) {
        if (border_nodes.contains(i)) {
            mirror_nodes[i] = draft.nodes[i];
        }
        else {
            mirror_nodes[i] = NodeInput::create(draft.nodes[i]->pos);
            mirror_nodes[i]->pos[axis] *= -1.0;
        }
    }

    draft.nodes.reserve(2 * n_nodes_prev);
    for (const auto& node: mirror_nodes) {
        add_node(node);
    }

    // Создать копии существующих ячеек с индексами новых вершин
    const id_t n_cells_prev = draft.cells.size();
    draft.cells.resize(2 * n_cells_prev);
    for (id_t i = n_cells_prev; i < draft.cells.size(); ++i) {
        draft.cells[i] = draft.cells[i - n_cells_prev];
        for (id_t& j: draft.cells[i].node_ids) {
            j = mirror_nodes[j]->m_builder.id;
        }
    }

    // Необходимо развернуть ячейки
    if (dimension() == 2) {
        // Разворачиваем полигоны, оставляя первую вершину на месте
        if (m_type == Type::TRI || m_type == Type::QUAD || m_type == Type::POLY) {
            for (id_t i = n_cells_prev; i < draft.cells.size(); ++i) {
                draft.cells[i].node_ids.push_back(draft.cells[i].node_ids.front());
                std::ranges::reverse(draft.cells[i].node_ids);
                draft.cells[i].node_ids.pop_back();
                std::ranges::reverse(draft.cells[i].face_bc);
            }
        }
        else {
            if (m_type != Type::AMR) {
                throw std::runtime_error("Grid::mirror: strange grid type");
            }

            throw std::runtime_error("Amr mirror not implemented");
        }
    }
    else {
        throw std::runtime_error("3D mirror not implemented");
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

    if (!m_draft) {
        throw std::runtime_error("Grid::rotate: no draft data");
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
        std::cerr << "Sorry, maybe later\n";
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
    DraftData& draft = *m_draft;

    int n_cells = draft.cells.size();
    draft.cells.resize(2 * n_cells);
    for (id_t i = 0; i < n_cells; ++i) {
        std::array node_ids{
            draft.cells[i].node_ids[0], draft.cells[i].node_ids[1],
            draft.cells[i].node_ids[2], draft.cells[i].node_ids[3]
        };

        std::array face_bc{
            draft.cells[i].face_bc[0], draft.cells[i].face_bc[1],
            draft.cells[i].face_bc[2], draft.cells[i].face_bc[3]
        };

        std::array vs = {
            draft.nodes[node_ids[0]]->pos, draft.nodes[node_ids[1]]->pos,
            draft.nodes[node_ids[2]]->pos, draft.nodes[node_ids[3]]->pos
        };

        id_t j = n_cells + i;
        draft.cells[i].type = CellType::TRIANGLE;
        draft.cells[j].type = CellType::TRIANGLE;

        if (choose_delaunay_diagonal(vs) == 0) {
            draft.cells[i].node_ids = {node_ids[0], node_ids[1], node_ids[2]};
            draft.cells[i].face_bc = {face_bc[0], face_bc[1], Boundary::INNER};

            draft.cells[j].node_ids = {node_ids[0], node_ids[2], node_ids[3]};
            draft.cells[j].face_bc = {Boundary::INNER, face_bc[2], face_bc[3]};
        }
        else {
            draft.cells[i].node_ids = {node_ids[0], node_ids[1], node_ids[3]};
            draft.cells[i].face_bc = {face_bc[0], Boundary::INNER, face_bc[3]};

            draft.cells[j].node_ids = {node_ids[1], node_ids[2], node_ids[3]};
            draft.cells[j].face_bc = {face_bc[1], face_bc[2], Boundary::INNER};
        }
    }
    m_type = Type::TRI;
}

void Grid::triangulation_quad_symmetry_() {
    DraftData& draft = *m_draft;

    // Создать центральные вершины
    std::vector<NodeInput::Ptr> central(draft.cells.size());
    for (id_t i = 0; i < draft.cells.size(); ++i) {
        Vector3d c = Vector3d::Zero();
        for (id_t j: draft.cells[i].node_ids) {
            c += draft.nodes[j]->pos;
        }
        c /= draft.cells[i].node_ids.size();
        central[i] = NodeInput::create(c);
    }
    for (auto node: central) {
        add_node(node);
    }

    int n_cells = draft.cells.size();
    draft.cells.resize(4 * n_cells);
    for (id_t i = 0; i < n_cells; ++i) {
        std::array node_ids{
            draft.cells[i].node_ids[0], draft.cells[i].node_ids[1],
            draft.cells[i].node_ids[2], draft.cells[i].node_ids[3]
        };

        std::array face_bc{
            draft.cells[i].face_bc[0], draft.cells[i].face_bc[1],
            draft.cells[i].face_bc[2], draft.cells[i].face_bc[3]
        };

        id_t central_id = central[i]->m_builder.id;

        id_t j1 = n_cells + i;
        id_t j2 = 2*n_cells + i;
        id_t j3 = 3*n_cells + i;

        draft.cells[i].type = CellType::TRIANGLE;
        draft.cells[i].node_ids = {node_ids[0], node_ids[1], central_id};
        draft.cells[i].face_bc = {face_bc[0], Boundary::INNER, Boundary::INNER};

        draft.cells[j1].type = CellType::TRIANGLE;
        draft.cells[j1].node_ids = {node_ids[1], node_ids[2], central_id};
        draft.cells[j1].face_bc = {face_bc[1], Boundary::INNER, Boundary::INNER};

        draft.cells[j2].type = CellType::TRIANGLE;
        draft.cells[j2].node_ids = {node_ids[2], node_ids[3], central_id};
        draft.cells[j2].face_bc = {face_bc[2], Boundary::INNER, Boundary::INNER};

        draft.cells[j3].type = CellType::TRIANGLE;
        draft.cells[j3].node_ids = {node_ids[3], node_ids[0], central_id};
        draft.cells[j3].face_bc = {face_bc[3], Boundary::INNER, Boundary::INNER};
    }
    m_type = Type::TRI;
}

void Grid::make_amr_2D_() {
    DraftData& draft = *m_draft;

    // Создать центральные вершины
    std::vector<NodeInput::Ptr> central(draft.cells.size());
    for (id_t i = 0; i < draft.cells.size(); ++i) {
        Vector3d c = Vector3d::Zero();
        for (id_t j: draft.cells[i].node_ids) {
            c += draft.nodes[j]->pos;
        }
        c /= draft.cells[i].node_ids.size();
        central[i] = NodeInput::create(c);
    }

    // Добавить центральные в массив
    for (const auto& node: central) {
        add_node(node);
    }

    struct face_t {
        id_t v1, v2;
        bool operator==(const face_t& other) const {
            return (v1 == other.v1 && v2 == other.v2) ||
                   (v1 == other.v2 && v2 == other.v1);
        }
    };
    struct face_hash {
        std::size_t operator()(const face_t& k) const noexcept {
            auto h1 = std::hash<id_t>()(k.v1);
            auto h2 = std::hash<id_t>()(k.v2);
            return (h1 + h2) * (h1 + h2 + 1) / 2 + std::min(h1, h2);
        }
    };

    // Уникальные грани
    std::unordered_map<const face_t, NodeInput::Ptr, face_hash> faces;
    for (auto& cell: draft.cells) {
        z_assert(cell.type == CellType::QUAD && cell.node_ids.size() == 4, "Bad cell");
        for (int i = 0; i < 4; ++i) {
            face_t face{cell.node_ids[i], cell.node_ids[(i + 1) % 4]};
            if (!faces.contains(face)) {
                Vector3d fc = 0.5*(draft.nodes[face.v1]->pos + draft.nodes[face.v2]->pos);
                faces[face] = NodeInput::create(fc);
            }
        }
    }

    // Добавить вершины граней в массив
    for (const auto& val: faces | std::views::values) {
        add_node(val);
    }

    // Преобразовать ячейки в AMR2D (заменяем полигоны на таблицы 3x3)
    for (id_t i = 0; i < draft.cells.size(); ++i) {
        auto& cell = draft.cells[i];
        face_t L{cell.node_ids[0], cell.node_ids[3]};
        face_t R{cell.node_ids[1], cell.node_ids[2]};
        face_t B{cell.node_ids[0], cell.node_ids[1]};
        face_t T{cell.node_ids[3], cell.node_ids[2]};

        z_assert(faces.contains(L), "Cell face not found");
        z_assert(faces.contains(R), "Cell face not found");
        z_assert(faces.contains(B), "Cell face not found");
        z_assert(faces.contains(T), "Cell face not found");

        cell.type = CellType::AMR2D;

        id_t vLB = cell.node_ids[0];
        id_t vRB = cell.node_ids[1];
        id_t vLT = cell.node_ids[3];
        id_t vRT = cell.node_ids[2];
        id_t v_C = central[i]->id();
        id_t v_L = faces[L]->id();
        id_t v_R = faces[R]->id();
        id_t v_B = faces[B]->id();
        id_t v_T = faces[T]->id();

        cell.node_ids = {
            vLB, v_B, vRB,
            v_L, v_C, v_R,
            vLT, v_T, vRT
        };

        // Переставить гран условия
        std::vector bc = draft.cells[i].face_bc;
        draft.cells[i].face_bc[Side2D::L] = bc[3];
        draft.cells[i].face_bc[Side2D::R] = bc[1];
        draft.cells[i].face_bc[Side2D::B] = bc[0];
        draft.cells[i].face_bc[Side2D::T] = bc[2];
    }
    m_type = Type::AMR;
}

void Grid::make_amr_3D_() {
    throw std::runtime_error("Grid::make_amr_3D_ not implemented");
}

void Grid::Grid::make_amr() {
    require_editable_();
    if (m_type != Type::QUAD) {
        throw std::runtime_error("Grid::make_amr: grid is not QUAD type");
    }
    if (!m_draft) {
        throw std::runtime_error("Grid::make_amr: no draft data");
    }
    if (dimension() == 2) {
        make_amr_2D_();
    }
    else {
        make_amr_3D_();
    }
}

void Grid::reduce_amr_2D_() {
    throw std::runtime_error("Grid::reduce_amr_2D_ not implemented");
}

void Grid::reduce_amr_3D_() {
    throw std::runtime_error("Grid::reduce_amr_3D_ not implemented");
}

void Grid::reduce_amr() {
    require_editable_();
    if (m_type != Type::AMR) {
        throw std::runtime_error("Grid::reduce_amr: grid is not AMR type");
    }
    if (!m_draft) {
        throw std::runtime_error("Grid::reduce_amr: no draft data");
    }
    if (dimension() == 2) {
        reduce_amr_2D_();
    }
    else {
        reduce_amr_3D_();
    }
}

} // namespace zephyr::geom