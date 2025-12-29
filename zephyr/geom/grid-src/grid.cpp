#include <limits>
#include <stdexcept>
#include <utility>

#include <zephyr/geom/grid.h>

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

struct CellArity {
    std::size_t n_nodes_min = 0;
    std::size_t n_nodes_max = 0; // 0 => "no upper bound"
    std::size_t n_faces_min = 0;
    std::size_t n_faces_max = 0; // 0 => "no upper bound"
    bool faces_equal_nodes = false; // for 2D polygon (and generally 2D cells)
};

CellArity arity_of(CellType t) {
    switch (t) {
        // 2D
        case CellType::TRIANGLE:   return {3, 3, 3, 3, true};
        case CellType::QUAD:       return {4, 4, 4, 4, true};
        case CellType::POLYGON:    return {3, 0, 3, 0, true};   // faces = edges = nodes

        case CellType::AMR2D:      return {9, 9, 4, 4, false};  // 3x3 nodes, faces like QUAD

            // 3D (VTK-like counts)
        case CellType::TETRA:      return {4, 4, 4, 4, false};
        case CellType::PYRAMID:    return {5, 5, 5, 5, false};
        case CellType::WEDGE:      return {6, 6, 5, 5, false};  // triangular prism: 5 faces
        case CellType::HEXAHEDRON: return {8, 8, 6, 6, false};

        case CellType::AMR3D:      return {27, 27, 6, 6, false}; // 3x3x3 nodes, faces like HEX

        case CellType::POLYHEDRON: return {0, 0, 0, 0, false};   // not supported here
        default:                   return {0, 0, 0, 0, false};
    }
}

void require_ptrs_non_null(const std::vector<NodeInput*>& nodes) {
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (!nodes[i]) {
            throw std::invalid_argument("add_cell: nodes contains null pointer at index " + std::to_string(i));
        }
    }
}

} // anonymous namespace

void Grid::BuildOptions::validate_or_throw(int dim) const {
    if (faces == FaceOption::none) {
        if (build_face_local_indices || build_twin_face || build_node_faces || compute_face_geometry) {
            throw std::invalid_argument("BuildOptions: face-related options require faces != none.");
        }
    }

    if (build_twin_face && faces != FaceOption::per_cell) {
        throw std::invalid_argument("BuildOptions: build_twin_face requires faces == per_cell.");
    }

    if (build_node_faces && faces == FaceOption::none) {
        throw std::invalid_argument("BuildOptions: build_node_faces requires faces != none.");
    }

    if (compute_face_geometry && faces == FaceOption::none) {
        throw std::invalid_argument("BuildOptions: compute_face_geometry requires faces != none.");
    }

    if (build_edges && dim < 3) {
        throw std::invalid_argument("BuildOptions: build_edges available only for 3D grids.");
    }
}

Grid::Grid() = default;

Grid::~Grid() = default;

Grid::Grid(Grid&& other) noexcept
    : m_state(other.m_state),
      m_draft(std::move(other.m_draft)),
      m_nodes(std::move(other.m_nodes)),
      m_cells(std::move(other.m_cells)),
      m_geom(std::move(other.m_geom)),
      m_faces(std::move(other.m_faces)),
      m_edges(std::move(other.m_edges)),
      m_inc(std::move(other.m_inc))
{
    // Ensure moved-from object stays in a safe, editable state.
    other.remove_options_();
}

Grid& Grid::operator=(Grid&& other) noexcept {
    if (this == &other) {
        return *this;
    }

    m_state = other.m_state;
    m_draft = std::move(other.m_draft);
    m_nodes = std::move(other.m_nodes);
    m_cells = std::move(other.m_cells);
    m_geom  = std::move(other.m_geom);
    m_faces = std::move(other.m_faces);
    m_edges = std::move(other.m_edges);
    m_inc   = std::move(other.m_inc);

    // Ensure moved-from object stays in a safe, editable state.
    other.remove_options_();

    return *this;
}

void Grid::remove_options_() {
    m_state = State::editable;
    m_draft.reset();
    m_geom.reset();
    m_faces.reset();
    m_edges.reset();
    m_inc.reset();
}

void Grid::DraftData::clear() {
    nodes.clear();
    cells.clear();
    node_buffer.clear();
    stamp = 1;
    next_node_id = 0;
}

inline std::uint32_t next_stamp(std::uint32_t s) noexcept {
    // 0 оставим "неиспользуемым" (на всякий случай), крутимся в [1..max]
    if (s == std::numeric_limits<std::uint32_t>::max()) {
        return 1u;
    }
    return s + 1u;
}

void Grid::require_editable_() const {
    if (m_state != Grid::State::editable) {
        throw std::logic_error("Grid is not editable. Call clear()/begin_build() first.");
    }
}
void Grid::require_finalized_() const {
    if (m_state != Grid::State::finalized) {
        throw std::logic_error("Grid is not finalized. Call finalize().");
    }
}

void Grid::clear() {
    // Возвращаемся в "чистое" editable состояние с нулевыми данными.
    m_state = State::editable;

    m_draft.reset();

    m_dim = 0;
    m_type = Type::NONE;
    m_options = {};
    m_nodes.clear();
    m_cells.clear();

    m_geom.reset();
    m_faces.reset();
    m_edges.reset();
    m_inc.reset();
}

void Grid::begin_build() {
    // Старт новой build-сессии: очищаем draft, но stamp инкрементируем,
    // чтобы старые NodeInput::m_builder.id/stamp считались "невалидными".
    m_state = State::editable;

    if (!m_draft) {
        m_draft.emplace(DraftData{});
    } else {
        // Очищаем буферы, сбрасываем id-счётчик
        m_draft->nodes.clear();
        m_draft->cells.clear();
        m_draft->next_node_id = 0;
        m_draft->stamp += 1;
    }

    // На всякий случай сбросим "финальную" часть
    m_dim = 0;
    m_type = Type::NONE;
    m_options = {};
    m_nodes.clear();
    m_cells.clear();

    m_geom.reset();
    m_faces.reset();
    m_edges.reset();
    m_inc.reset();
}

void Grid::reserve_nodes(id_t n_nodes) {
    if (!m_draft) {
        begin_build();
    }
    m_draft->nodes.reserve(n_nodes);
}

id_t Grid::add_node(NodeInput* v) {
    require_editable_();

    if (!v) {
        throw std::invalid_argument("add_node: NodeInput pointer is null.");
    }
    if (!m_draft) {
        begin_build();
    }

    auto& bs = v->m_builder;
    if (bs.stamp == m_draft->stamp && bs.id != invalid_id) {
        return bs.id;
    }

    // Первый раз видим этот указатель в текущей build-сессии
    bs.id = m_draft->next_node_id++;
    bs.stamp = m_draft->stamp;

    z_assert(m_draft->nodes.size() == bs.id, "Grid::add_node(): nodes must be unique.");

    m_draft->nodes.push_back(v); // index == id
    return bs.id;
}

void Grid::reserve_cells(id_t n_cells) {
    if (!m_draft) {
        begin_build();
    }
    m_draft->cells.reserve(n_cells);
}

id_t Grid::add_cell(CellType type,
                    const std::vector<NodeInput*>& nodes,
                    const std::vector<Boundary>& face_bc) {
    require_editable_();

    // POLYHEDRON is explicitly not supported by this overload
    if (type == CellType::POLYHEDRON) {
        throw std::invalid_argument(
            "add_cell(POLYHEDRON): not supported. Use add_polyhedron(nodes, faces, faces_bc).");
    }

    // Ensure draft exists and set/check dimension.
    if (!m_draft.has_value()) {
        begin_build(); // must create m_draft
        m_dim = get_dimension(type);
        m_type = grid_type(type);
    } else {
        const int dim = get_dimension(type);
        if (m_dim == 0) {
            m_dim = dim;
        } else if (m_dim != dim) {
            throw std::runtime_error("add_cell: wrong dimension (cell dim mismatch with grid dim).");
        }
        m_type = promotion(m_type, type);
    }

    // Validate arity (nodes count and face_bc count).
    const CellArity a = arity_of(type);
    if (a.n_nodes_min == 0) {
        throw std::invalid_argument("add_cell: unsupported CellType.");
    }

    const std::size_t n_nodes = nodes.size();
    if (n_nodes < a.n_nodes_min || (a.n_nodes_max != 0 && n_nodes != a.n_nodes_max)) {
        throw std::invalid_argument("add_cell: invalid number of nodes for this CellType.");
    }

    const std::size_t n_faces = a.faces_equal_nodes ? n_nodes : a.n_faces_min;

    if (!face_bc.empty() && face_bc.size() != n_faces) {
        throw std::invalid_argument("add_cell: face_bc must be empty or have correct size for this CellType.");
    }

    require_ptrs_non_null(nodes);

    // Build DraftCell
    DraftCell dc{.type = type};
    dc.node_ids.reserve(n_nodes);

    for (NodeInput* p : nodes) {
        dc.node_ids.push_back(add_node(p)); // pointer-identity dedup
    }

    // BCs: if not provided, fill with UNDEFINED
    if (face_bc.empty()) {
        dc.face_bc.assign(n_faces, Boundary::UNDEFINED);
    } else {
        dc.face_bc = face_bc;
    }

    auto& d = *m_draft;
    const id_t cell_id = static_cast<id_t>(d.cells.size());
    d.cells.push_back(std::move(dc));
    return cell_id;
}

id_t Grid::add_polyhedron(const std::vector<NodeInput*>& nodes,
                          const std::vector<std::vector<NodeInput*>>& faces,
                          const std::vector<Boundary>& faces_bc) {
    require_editable_();

    // draft init + dim set/check
    if (!m_draft) {
        begin_build();
        m_dim = 3;
        m_type = Type::POLYHEDRON;
    } else {
        if (m_dim != 3) {
            throw std::runtime_error("add_polyhedron: wrong dimension (grid dim mismatch).");
        }
        if (m_type != Type::POLYHEDRON) {
            throw std::runtime_error("add_cell: cell type mismatch with grid type.");
        }
    }

    if (nodes.empty()) {
        throw std::invalid_argument("add_polyhedron: nodes is empty.");
    }
    if (faces.empty()) {
        throw std::invalid_argument("add_polyhedron: faces is empty.");
    }
    if (!faces_bc.empty() && faces_bc.size() != faces.size()) {
        throw std::invalid_argument("add_polyhedron: faces_bc must be empty or size==faces.size().");
    }

    // Build DraftCell
    DraftCell dc{.type=CellType::POLYHEDRON};

    // Convert cell nodes to global ids
    dc.node_ids.reserve(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        NodeInput* p = nodes[i];
        if (!p) {
            throw std::invalid_argument("add_polyhedron: nodes contains null pointer at index " + std::to_string(i));
        }
        dc.node_ids.push_back(add_node(p));
    }

    // Validate duplicates inside polyhedron node list (almost always a bug)
    {
        auto tmp = dc.node_ids;
        std::sort(tmp.begin(), tmp.end());
        if (std::adjacent_find(tmp.begin(), tmp.end()) != tmp.end()) {
            throw std::invalid_argument("add_polyhedron: nodes contain duplicates (by pointer/global id).");
        }
    }

    // Fill face BCs
    if (faces_bc.empty()) {
        dc.face_bc.assign(faces.size(), Boundary::UNDEFINED);
    } else {
        dc.face_bc = faces_bc;
    }

    // Build polyhedron face CSR: offsets + flattened face nodes (global ids)
    dc.poly_face_offsets.resize(faces.size() + 1);
    dc.poly_face_offsets[0] = 0;

    // For membership check: cell node ids sorted
    std::vector<id_t> cell_nodes_sorted = dc.node_ids;
    std::sort(cell_nodes_sorted.begin(), cell_nodes_sorted.end());

    auto cell_contains = [&](id_t nid) -> bool {
        return std::binary_search(cell_nodes_sorted.begin(), cell_nodes_sorted.end(), nid);
    };

    std::size_t total_face_nodes = 0;
    for (const auto& f : faces) total_face_nodes += f.size();
    dc.poly_face_nodes.reserve(total_face_nodes);

    id_t cursor = 0;
    for (std::size_t fi = 0; fi < faces.size(); ++fi) {
        const auto& fnodes = faces[fi];
        if (fnodes.size() < 3) {
            throw std::invalid_argument("add_polyhedron: face " + std::to_string(fi) + " has < 3 nodes.");
        }

        // Convert face nodes to global ids
        std::vector<id_t> face_ids;
        face_ids.reserve(fnodes.size());

        for (std::size_t k = 0; k < fnodes.size(); ++k) {
            NodeInput* p = fnodes[k];
            if (!p) {
                throw std::invalid_argument("add_polyhedron: face has null node pointer at face " +
                                            std::to_string(fi) + ", k=" + std::to_string(k));
            }
            const id_t gid = add_node(p);

            // Enforce that face nodes belong to the polyhedron node set.
            if (!cell_contains(gid)) {
                throw std::invalid_argument(
                    "add_polyhedron: face node is not listed in 'nodes' (face " + std::to_string(fi) + ").");
            }

            face_ids.push_back(gid);
        }

        // Detect duplicates inside a face (usually invalid)
        {
            auto tmp = face_ids;
            std::sort(tmp.begin(), tmp.end());
            if (std::adjacent_find(tmp.begin(), tmp.end()) != tmp.end()) {
                throw std::invalid_argument("add_polyhedron: face " + std::to_string(fi) + " has duplicate nodes.");
            }
        }

        // Append to flattened storage
        for (id_t gid : face_ids) {
            dc.poly_face_nodes.push_back(gid);
            ++cursor;
        }
        dc.poly_face_offsets[fi + 1] = cursor;
    }

    // Optional: basic sanity for polyhedron (not strict, but useful)
    if (dc.node_ids.size() < 4) {
        throw std::invalid_argument("add_polyhedron: polyhedron must have at least 4 unique nodes.");
    }
    if (faces.size() < 4) {
        throw std::invalid_argument("add_polyhedron: polyhedron must have at least 4 faces.");
    }

    // Store draft cell
    auto& d = *m_draft;
    const id_t cell_id = static_cast<id_t>(d.cells.size());
    d.cells.push_back(std::move(dc));
    return cell_id;
}

void GeometryCache::clear() {
    face_centroids.clear();
    face_normals.clear();
    face_areas.clear();
    cell_centroids.clear();
    cell_volumes.clear();
}

void FacesCache::clear() {
    cell_faces.clear();
    face_nodes.clear();
    owner_cell.clear();
    neighbor_cell.clear();
    face_bc.clear();
    twin_face.clear();
    local_nodes.clear();
}

std::size_t FacesCache::n_faces() const noexcept {
    return owner_cell.size();
}

/// @brief Basic consistency check (throws on mismatch).
void FacesCache::validate_or_throw(std::size_t n_cells_expected) const {
    if (owner_cell.size() != neighbor_cell.size() || owner_cell.size() != face_bc.size()) {
        throw std::logic_error("FacesCache: owner/neighbor/bc arrays size mismatch.");
    }
    if (cell_faces.size() != n_cells_expected) {
        throw std::logic_error("FacesCache: cell_faces row count mismatch.");
    }
    if (face_nodes.size() != owner_cell.size()) {
        throw std::logic_error("FacesCache: face_nodes row count mismatch.");
    }
    if (!twin_face.empty() && twin_face.size() != owner_cell.size()) {
        throw std::logic_error("FacesCache: twin_face size mismatch.");
    }
    if (!local_nodes.offsets.empty() && local_nodes.size() != owner_cell.size()) {
        throw std::logic_error("FacesCache: local_nodes row count mismatch.");
    }
}

void EdgesCache::clear() {
    edges.clear();
    edge_faces.clear();
    edge_cells.clear();
}

std::size_t EdgesCache::n_edges() const noexcept {
    return edges.size();
}

void EdgesCache::validate_or_throw() const {
    if (edges.size() != edge_faces.size() || edges.size() != edge_cells.size()) {
        throw std::logic_error("EdgesCache: arrays size mismatch.");
    }
}

void IncidenceCache::clear() {
    node_cells.clear();
    node_faces.clear();
}

bool Grid::faces_per_cell() const noexcept {
    return m_options.faces == BuildOptions::FaceOption::per_cell;
}

bool Grid::has_faces() const noexcept {
    return m_faces.has_value();
}

std::size_t Grid::n_faces() const noexcept {
    return has_faces() ? m_faces->n_faces() : 0;
}

std::size_t Grid::total_nodes_per_cell() const noexcept {
    size_t count = 0;
    for (auto cell: m_cells) {
        count += cell.nodes.size();
    }
    return count;
}

bool Grid::has_cells_geometry() const noexcept {
    return m_geom.has_value() && !m_geom->cell_centroids.empty();
}

bool Grid::has_faces_geometry() const noexcept {
    return m_geom.has_value() && !m_geom->face_centroids.empty();
}

} // namespace zephyr::geom