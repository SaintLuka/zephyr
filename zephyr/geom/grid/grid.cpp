#include <limits>
#include <stdexcept>
#include <utility>
#include <ranges>

#include <zephyr/geom/grid.h>
#include <zephyr/geom/indexing.h>

#include <zephyr/geom/primitives/line.h>
#include <zephyr/geom/primitives/triangle.h>
#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/geom/primitives/quad.h>
#include <zephyr/geom/primitives/cube.h>
#include <zephyr/geom/primitives/polyhedron.h>

#include <zephyr/utils/threads.h>

using zephyr::utils::threads;

namespace zephyr::geom {

namespace {

// Grid type by cell type
Grid::Type grid_type(CellType type) {
    if (type == CellType::AMR2D || type == CellType::AMR3D) {
        return Grid::Type::AMR;
    }
    if (type == CellType::TRIANGLE || type == CellType::TETRA) {
        return Grid::Type::TRI;
    }
    if (type == CellType::QUAD || type == CellType::HEXAHEDRON) {
        return Grid::Type::QUAD;
    }
    if (type == CellType::POLYHEDRON) {
        return Grid::Type::POLYHEDRON;
    }
    // Polygons and standard 3d elements (Wedge, Pyramid)
    return Grid::Type::POLY;;
}

// New grid type by appended cell (doesn't check dimensions)
Grid::Type promotion(Grid::Type type, CellType cell) {
    // First time, always promote
    if (type == Grid::Type::NONE) {
        return grid_type(cell);
    }

    // AMR type doesn't change
    if (type == Grid::Type::AMR) {
        if (cell == CellType::AMR2D || cell == CellType::AMR3D) {
            return type;
        }
        throw std::runtime_error("bad promotion #1: inconsistent cell type and grid type.");
    }
    // type != Grid::Type::AMR
    if (cell == CellType::AMR2D || cell == CellType::AMR3D) {
        throw std::runtime_error("bad promotion #2: inconsistent cell type and grid type.");
    }

    // General POLYHEDRON always promote to Type::HARD
    if (cell == CellType::POLYHEDRON) {
        return Grid::Type::POLYHEDRON;
    }

    // Can promote Type::TRI to Type::POLY or Type::HARD
    if (type == Grid::Type::TRI) {
        return (cell == CellType::TRIANGLE || cell == CellType::TETRA) ? type : Grid::Type::POLY;
    }

    // Can promote Type::QUAD to Type::POLY or Type::HARD
    if (type == Grid::Type::QUAD) {
        return (cell == CellType::QUAD || cell == CellType::HEXAHEDRON) ? type : Grid::Type::POLY;
    }

    return Grid::Type::POLY;
}

} // anonymous namespace

Grid::Grid() = default;

Grid::~Grid() = default;

Grid::Grid(Grid&& other) noexcept
    : m_state(other.m_state),
      m_nodes(std::move(other.m_nodes)),
      m_cells(std::move(other.m_cells))
{
    // Ensure moved-from object stays in a safe, editable state.
    other.clear();
}

Grid& Grid::operator=(Grid&& other) noexcept {
    if (this == &other) {
        return *this;
    }

    m_state = other.m_state;
    m_dim = other.m_dim;
    m_type = other.m_type;
    m_nodes = std::move(other.m_nodes);
    m_cells = std::move(other.m_cells);

    // Ensure moved-from object stays in a safe, editable state.
    other.clear();

    return *this;
}

void Grid::require_editable_() const {
    if (m_state != State::Editable) {
        throw std::logic_error("Grid is not editable. Call clear() first.");
    }
}
void Grid::require_finalized_() const {
    if (m_state != State::Finalized) {
        throw std::logic_error("Grid is not finalized. Call finalize().");
    }
}

void Grid::clear() {
    // Возвращаемся в "чистое" editable состояние с нулевыми данными.
    m_state = State::Editable;
    m_dim = 0;
    m_type = Type::NONE;
    m_nodes.clear();
    m_cells.clear();
}

void Grid::reserve_nodes(id_t n_nodes) {
    m_nodes.reserve(n_nodes);
}

id_t Grid::add_node(Node::Ref node) {
    require_editable_();

    if (!node) {
        throw std::invalid_argument("add_node: Node pointer is null.");
    }

    if (node->id() != invalid_id) {
        if (m_nodes.size() < node->id()) {
            throw std::runtime_error("add_node: existed Node index is too large?");
        }
        if (m_nodes[node->id()] != *node) {
            throw std::runtime_error("add_node: another Node already exists.");
        }
        return node->id();
    }

    // Первый раз видим этот указатель
    node->m_id = m_nodes.size();
    m_nodes.push_back(*node);
    return node->id();
}

void Grid::reserve_cells(id_t n_cells) {
    m_cells.reserve(n_cells);
}

std::vector<id_t> Grid::add_nodes(const std::vector<Node::Ptr>& nodes) {
    std::vector<id_t> node_ids(nodes.size());
    for (int i = 0; i < nodes.size(); ++i) {
        if (!nodes[i]) {
            throw std::invalid_argument("add_cell: node is nullptr at index " + std::to_string(i));
        }
        node_ids[i] = add_node(nodes[i]);
    }
    node_ids.reserve(nodes.size());

    // Check duplicates
    auto tmp = node_ids;
    std::ranges::sort(tmp);
    if (std::ranges::adjacent_find(tmp) != tmp.end()) {
        throw std::invalid_argument("add_polyhedron: nodes contain duplicates (by pointer/global id).");
    }

    return node_ids;
}

id_t Grid::add_cell(CellType type,
                    const std::vector<Node::Ptr>& nodes,
                    const std::vector<Boundary>& faces_bc) {
    require_editable_();

    // POLYHEDRON is explicitly not supported by this overload
    if (type == CellType::POLYHEDRON) {
        throw std::invalid_argument(
            "add_cell(POLYHEDRON): not supported. Use add_polyhedron(nodes, faces, faces_bc).");
    }

    // Ensure draft exists and set/check dimension.
    const int dim = get_dimension(type);
    if (m_dim == 0) {
        m_dim = dim;
    } else if (m_dim != dim) {
        throw std::runtime_error("add_cell: wrong dimension (cell dim mismatch with grid dim).");
    }
    m_type = promotion(m_type, type);

    // Добавить узлы, если их нет
    auto node_ids = add_nodes(nodes);

    // Добавить ячейку
    const id_t cell_id = m_cells.size();
    m_cells.emplace_back(type, std::move(node_ids));
    if (!faces_bc.empty()) {
        // Проинициализирует грани
        m_cells.back().set_face_bc(faces_bc);
    }
    return cell_id;
}

id_t Grid::add_polyhedron(const std::vector<Node::Ptr>& nodes,
                          const std::vector<std::vector<Node::Ptr>>& faces,
                          const std::vector<Boundary>& faces_bc) {
    require_editable_();

    if (m_cells.empty()) {
        m_dim = 3;
        m_type = Type::POLYHEDRON;
    } else {
        if (m_dim != 3) {
            throw std::runtime_error("Grid::add_polyhedron: wrong dimension (grid dim mismatch).");
        }
        m_type = Type::POLYHEDRON;
    }

    if (nodes.empty()) {
        throw std::invalid_argument("Grid::add_polyhedron: nodes is empty.");
    }
    if (faces.empty()) {
        throw std::invalid_argument("Grid::add_polyhedron: faces is empty.");
    }
    if (!faces_bc.empty() && faces_bc.size() != faces.size()) {
        throw std::invalid_argument("Grid::add_polyhedron: faces_bc must be empty or size==faces.size().");
    }

    // Добавить узлы, если их нет
    auto node_ids = add_nodes(nodes);

    // Добавить ячейку
    const id_t cell_id = m_cells.size();
    m_cells.emplace_back(CellType::POLYHEDRON, std::move(node_ids));

    // Надо собрать локальные индексы
    std::unordered_map<id_t, int> glob2loc;
    for (int i = 0; i < nodes.size(); ++i) {
        glob2loc[nodes[i]->id()] = i;
    }

    std::vector<std::vector<int>> local_indices(faces.size());
    for (int i = 0; i < faces.size(); ++i) {
        local_indices[i].reserve(faces[i].size());
        for (int j = 0; j < faces[i].size(); ++j) {
            local_indices[i].push_back(glob2loc[faces[i][j]->id()]);
        }
    }

    m_cells.back().set_faces(local_indices);
    if (!faces_bc.empty()) {
        m_cells.back().set_face_bc(faces_bc);
    }
    return cell_id;
}

bool Grid::has_faces() const noexcept {
    return !m_cells.empty() && m_cells[0].has_faces();
}

std::size_t Grid::total_faces_per_cell() const noexcept {
    if (m_type == Type::TRI) {
        return (m_dim == 2 ? 3 : 4) * n_cells();
    }
    if (m_type == Type::QUAD || m_type == Type::AMR) {
        return (m_dim == 2 ? 4 : 6) * n_cells();
    }
    size_t counter = 0;
    for (const auto& cell: m_cells) {
        counter += cell.n_faces();
    }
    return counter;
}

std::size_t Grid::total_nodes_per_cell() const noexcept {
    if (m_type == Type::TRI) {
        return (m_dim == 2 ? 3 : 4) * n_cells();
    }
    if (m_type == Type::QUAD) {
        return (m_dim == 2 ? 4 : 8) * n_cells();
    }
    if (m_type == Type::AMR) {
        return (m_dim == 2 ? 9 : 27) * n_cells();
    }

    size_t count = 0;
    for (const auto& cell: m_cells) {
        count += cell.n_nodes();
    }
    return count;
}

void Grid::initialize_faces_(const BuildOptions& options) {
    auto init_faces = [](Cell& cell) {
        if (cell.has_faces()) return;
        if (cell.type() == CellType::POLYHEDRON) {
            throw std::runtime_error("Grid::initialize_faces_: POLYHEDRON without faces");
        }
        cell.init_faces();
    };

    threads::for_each(m_cells.begin(), m_cells.end(), init_faces);
}

void Grid::initialize_geom_(const BuildOptions& options) {
    threads::for_each(m_cells.begin(), m_cells.end(),
        [this](Cell& cell) {
            cell.calc_geom(m_nodes);
        });
}

void Grid::finalize(const BuildOptions& options) {
    require_editable_();

    if (m_cells.empty() || m_nodes.empty()) {
        throw std::runtime_error("Grid::finalize: no draft data");
    }

    // Инициализация граней
    if (options.build_faces) {
        initialize_faces_(options);
    }

    // Инициализация геометрии
    initialize_geom_(options);

    // Найти соседей через грань
    if (options.build_faces) {
        // Сделаем "лениво" через хэш граней (face hash, adjacent (cell, face_id) )
        std::unordered_map<FaceKey, std::vector<std::pair<id_t, int>>> faces;
        for (id_t ic = 0; ic < m_cells.size(); ++ic) {
            const auto& cell = m_cells[ic];
            if (!cell.has_faces()) {
                throw std::runtime_error("Grid::finalize: cell without faces");
            }
            for (int iface = 0; iface < cell.n_faces(); ++iface) {
                const Face& face = cell.get_face(iface);
                FaceKey face_key{cell, face};
                faces[face_key].emplace_back(ic, iface);
            }
        }

        // Проставить соседнюю ячейку
        for (const auto& adjacent: faces | std::views::values) {
            if (adjacent.empty()) {
                throw std::runtime_error("Grid::finalize: face without adjacent cells");
            }
            if (adjacent.size() == 1) {
                // Граница
                auto [ic, iface] = adjacent.front();
                if (m_cells[ic].get_face(iface).bc() == Boundary::INNER) {
                    m_cells[ic].set_bc(iface, Boundary::UNDEFINED);
                }
                continue;
            }
            if (adjacent.size() == 2) {
                auto [ic, iface] = adjacent[0];
                auto [jc, jface] = adjacent[1];
                m_cells[ic].set_neib(iface, jc, jface);
                m_cells[jc].set_neib(jface, ic, iface);
                continue;
            }
            throw std::runtime_error("Grid::finalize: more than two cells link to one face");
        }
    }

    m_state = State::Finalized;
}

} // namespace zephyr::geom