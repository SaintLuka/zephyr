#pragma once

#include <memory>
#include <vector>
#include <optional>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/boundary.h>
#include <zephyr/geom/cell_type.h>
#include <zephyr/utils/csr.h>

namespace zephyr::geom {

// ============================================================================
//                          Basic types & utilities
// ============================================================================

using id_t = std::uint32_t;
static constexpr id_t invalid_id = static_cast<id_t>(-1);

/// @brief Draft node object used only during grid construction.
/// User creates/moves/owns these objects and passes pointers to Grid.
/// Grid uses hidden builder state for pointer-identity deduplication.
struct GridNode {
    using Ptr = std::shared_ptr<GridNode>;
    using Ref = const std::shared_ptr<GridNode>&;

    Vector3d pos{0.0, 0.0, 0.0};
    Boundary bc{Boundary::UNDEFINED};

    static Ptr create(const Vector3d& v) {
        return std::make_shared<GridNode>(v);
    }

    static Ptr create(double x, double y) {
        return std::make_shared<GridNode>(Vector3d{x, y, 0.0});
    }

    GridNode() = default;

    GridNode(const Vector3d& v) : pos(v) { }

    id_t id() const noexcept { return m_id; }

private:
    friend class Grid;

    // Builder-only state
    id_t m_id{invalid_id};
};

// ============================================================================
//               Finalized entities, minimal Grid description
// ============================================================================

/// @brief Finalized unique node stored in Grid after finalize().
struct Node {
    Vector3d pos{0.0, 0.0, 0.0};
    Boundary bc{Boundary::UNDEFINED};
};

/// @brief Finalized cell stored in Grid after finalize().
struct Cell {
    CellType type{CellType::QUAD};
    std::vector<id_t> nodes{}; ///< node ids (ordering matters, VTK)
};

// ============================================================================
//                  Geometry and topology caches (optional)
// ============================================================================

/// @brief Optional computed geometry for faces/cells.
/// Presence depends on finalize options.
struct GeometryCache {
    // Cell geometry (size = n_cells)
    std::vector<Vector3d> cell_centroids{};
    std::vector<double>   cell_volumes{};

    // Face geometry (size = n_faces)
    std::vector<Vector3d> face_centroids{};
    std::vector<Vector3d> face_normals{};   // convention: outward from owner cell
    std::vector<double>   face_areas{};

    /// @brief Clear all cached geometry.
    void clear();
};

/// @brief Faces-derived data (optional). Built in finalize() depending on BuildOptions.
/// Invariants (when built):
///   - n_faces = owner_cell.size() = neighbor_cell.size() = face_bc.size()
///   - face_nodes.size() == n_faces
///   - cell_faces.size() == n_cells
struct FacesCache {
    /// @brief For each cell: list of face ids (size = n_cells rows).
    Csr<id_t, id_t> cell_faces{};

    /// @brief For each face: owner cell id (size = n_faces).
    std::vector<id_t> owner_cell{};

    /// @brief For each face: neighbor cell id (size = n_faces).
    /// invalid_id for boundary faces.
    std::vector<id_t> neighbor_cell{};

    /// @brief For each face: boundary tag/type (size = n_faces).
    /// Meaningful for boundary faces (neighbor_cell == invalid_id).
    std::vector<Boundary> face_bc{};

    /// @brief For each face: list of global node ids (size = n_faces).
    Csr<id_t, id_t> face_nodes{};

    // ---- Optional connectivity extras ----

    /// @brief Opposite oriented face id for per-cell storage (size = n_faces).
    /// invalid_id for boundary faces or if not built / not applicable.
    std::vector<id_t> twin_face{};

    /// @brief Local node indices in owner cell for each face (size = n_faces).
    Csr<std::uint16_t, id_t> local_nodes{};

    /// @brief Clear all face caches.
    void clear();

    /// @brief Number of faces in cache.
    std::size_t n_faces() const noexcept;

    /// @brief Basic consistency check (throws on mismatch).
    void validate_or_throw(std::size_t n_cells_expected) const;
};

/// @brief Edges-derived data (optional). Built in finalize() depending
/// on BuildOptions, size = n_edges.
struct EdgesCache {
    /// @brief Pairs of nodes global indices
    std::vector<std::array<id_t, 2>> edges{};

    /// @brief Optional computed incidence for edges
    Csr<id_t, id_t> edge_faces{};
    Csr<id_t, id_t> edge_cells{};

    /// @brief Clear all face caches.
    void clear();

    /// @brief Number of faces in cache.
    std::size_t n_edges() const noexcept;

    /// @brief Basic consistency check (throws on mismatch).
    void validate_or_throw() const;
};

/// @brief Optional computed incidence for nodes.
struct IncidenceCache {
    Csr<id_t, id_t> node_cells;
    Csr<id_t, id_t> node_faces;

    void clear();
};

// ============================================================================
//                                  Grid
// ============================================================================

/// @brief Intermediate Grid representation for mesh generation and conversion.
///
/// Grid exists in two states:
///
/// 1) Editable (minimal) state:
///    - Created by adding nodes (as pointers to NodeInput), and then adding 
///      cells as lists of NodeInput pointers.
///    - Node deduplication is by pointer identity (same NodeInput object => same index),
///      implemented in O(1) by storing a hidden builder id/stamp inside Node.
///    - In this state, Grid holds only minimal data needed to describe topology:
///      unique nodes (as references to inputs) and cell-to-node connectivity.
///
/// 2) Finalized (extended) state:
///    - Produced by calling finalize(options).
///    - Depending on options, Grid builds faces, adjacency, optional node incidence 
///      (node->cells/faces), and optional geometry caches (areas, volumes, etc).
///    - After finalize(), the editable draft storage is released to save memory,
///      and the Grid becomes read-only with respect to topology.
///
/// Typical pipeline:
///   Generator -> Grid (editable) -> grid.finalize(build_options)
///             -> solver-specific mesh conversion.
class Grid {
public:
    /// @brief Grid state: editable (draft) or finalized (extended).
    enum class State {
        editable, finalized
    };

    /// @brief Grid type
    enum class Type {
        NONE,  ///< Undefined
        AMR,   ///< Adaptive 2D or 3D cells
        TRI,   ///< Only Triangle (2D) or Tetra (3D) cells
        QUAD,  ///< Only Quad (2D) or Hex (3D) cells
        POLY,  ///< Common polygons (2D) or classic 3D elements
        POLYHEDRON   ///< Common polyhedra (3D)
    };

    /// @brief Options controlling what finalize() computes/builds.
    struct BuildOptions {
        /// @brief Controls face construction and storage mode.
        enum class FaceOption {
            none,     ///< Do not build faces cache at all.
            unique,   ///< Build one face per topological face (internal faces are shared).
            per_cell  ///< Build faces per cell (internal faces duplicated, oriented per owner).
        };

        FaceOption faces{FaceOption::none};

        bool build_face_local_indices{false}; ///< Build FacesCache::local_nodes (per-face local indices in owner cell).
        bool build_twin_face{false};          ///< Build FacesCache::twin_face (only meaningful for per_cell).

        bool build_edges{false};              ///< Build EdgesCache

        bool build_node_cells{false};         ///< Build IncidenceCache::node_cells (CSR).
        bool build_node_faces{false};         ///< Build IncidenceCache::node_faces (CSR). Requires faces != none.

        bool compute_face_geometry{false};    ///< Build GeometryCache face arrays. Requires faces != none.
        bool compute_cell_geometry{false};    ///< Build GeometryCache cell arrays.

        /// @brief Validate option compatibility.
        void validate_or_throw(int dim) const;
    };

public:
    Grid();
    ~Grid();

    Grid(const Grid&) = delete;
    Grid& operator=(const Grid&) = delete;
    Grid(Grid&&) noexcept;
    Grid& operator=(Grid&&) noexcept;

    /// @brief Access current state.
    State state() const noexcept { return m_state; }

    /// @brief True if grid is in editable state.
    bool is_editable() const noexcept { return m_state == State::editable; }

    /// @brief True if grid is finalized.
    bool is_finalized() const noexcept { return m_state == State::finalized; }

    // --------------------------
    // Editable API
    // --------------------------

    /// @brief Reset grid to empty editable state (invalidates all data).
    void clear();

    /// @brief Begin a new editable build session (increments internal stamp).
    ///        Useful when reusing the same Vertex objects across builds.
    void begin_build();

    /// @brief Reserve array for nodes
    void reserve_nodes(id_t n_nodes);

    /// @brief Reserve array for cells
    void reserve_cells(id_t n_cells);

    /// @brief Add a node via shared pointer; returns node index.
    id_t add_node(GridNode::Ref node);

    /// @brief Add a cell using a vector of
    id_t add_cell(CellType type, const std::vector<GridNode::Ptr>& nodes,
                  const std::vector<Boundary>& face_bc = {});

    /// @brief Add a cell of special kind (CellType::POLYHEDRON)
    id_t add_polyhedron(const std::vector<GridNode::Ptr>& nodes,
                        const std::vector<std::vector<GridNode::Ptr>>& faces,
                        const std::vector<Boundary>& faces_bc = {});

    // --------------------------
    // Grid Modification
    // --------------------------

    /// @brief Move grid to vector shift
    void move(const Vector3d& shift);

    /// @brief Scale grid
    void scale(double q);

    /// @brief Rotate grid by angle
    void rotate(double phi);

    /// @brief Rotate grid by rotation matrix R
    void rotate(const Matrix3d& R);

    /// @brief Преобразовать узлы сетки
    void transform(std::function<Vector3d(const Vector3d&)>& func);

    /// @brief Достроить зеркальное отражение сетки. Все узлы сетки должны
    /// лежать в полуплоскости X < 0 или X > 0.
    void mirror_x();

    /// @brief Достроить зеркальное отражение сетки. Все узлы сетки должны
    /// лежать в полуплоскости Y < 0 или Y > 0.
    void mirror_y();

    /// @brief Достроить зеркальное отражение сетки. Все узлы сетки должны
    /// лежать в полуплоскости Z < 0 или Z > 0.
    void mirror_z();

    /// @brief Выполнить триангуляцию существующей QUAD сетки.
    /// @param mode 1: триангуляция каждого квадрата на две части,
    ///             2: триангуляция каждого квадрата на 4 части
    ///                (с добавлением центральной точки)
    void triangulation(int mode);

    /// @brief Делает сетку из кубиков сеткой из пирамидок
    void pyramidize();

    /// @brief Convert Type::QUAD grid to Type::AMR grid (add additional unique nodes)
    void make_amr();

    /// @brief Convert Type::AMR grid to Type::QUAD grid (remove intermediate nodes)
    void reduce_amr();

    // --------------------------
    // Finalize API
    // --------------------------

    /// @brief Build requested extended connectivity/geometry and switch to finalized state.
    ///        Releases draft storage unless keep_input_nodes == true.
    void finalize(const BuildOptions& opt);

    // --------------------------
    // Read-only access (finalized)
    // --------------------------

    /// @brief Grid dimension (2 or 3)
    int dimension() const noexcept { return m_dim; }

    /// @brief Adaptive mesh (AMR cells is used)
    bool adaptive() const noexcept { return m_type == Type::AMR; }

    /// @brief 3D polyhedral mesh?
    bool polyhedral() const { throw std::runtime_error("sdfsdf"); }

    /// @brief Number of unique nodes (finalized).
    std::size_t n_nodes() const noexcept { return m_nodes.size(); }

    /// @brief Number of cells (finalized).
    std::size_t n_cells() const noexcept { return m_cells.size(); }

    /// @brief Access finalized nodes.
    const std::vector<Node>& nodes() const noexcept { return m_nodes; }

    /// @brief Access finalized cells.
    const std::vector<Cell>& cells() const noexcept { return m_cells; }

    /// @brief Faces were built  per cell (internal
    /// faces duplicated, oriented per owner).
    bool faces_per_cell() const noexcept;

    /// @brief Faces is computed?
    bool has_faces() const noexcept;

    /// @brief Number of faces if computed (else 0)
    std::size_t n_faces() const noexcept;

    const FacesCache& faces() const noexcept { return *m_faces; }

    /// @brief Number of nodes for each cell (with duplicates).
    std::size_t total_nodes_per_cell() const noexcept;

    bool has_geometry() const noexcept { return m_geom.has_value(); }
    bool has_cells_geometry() const noexcept;
    bool has_faces_geometry() const noexcept;

    const GeometryCache& geometry() const noexcept { return *m_geom; }


    // --------------------------
    // Validation / utilities
    // --------------------------

    /// @brief Validate basic invariants (ids in range, adjacency consistency, etc.).
    /// @return true if valid, false otherwise; optionally fills report.
    bool validate(std::string* report) const;

    bool validate_draft(std::string* report = nullptr) const;
    bool validate_finalized_basic(std::string* report = nullptr) const;
    bool validate_finalized_full(std::string* report = nullptr) const;

private:
    // ==========================
    // Draft (editable) storage
    // ==========================

    struct DraftCell {
        CellType type{CellType::QUAD};
        std::vector<id_t> node_ids{};
        std::vector<Boundary> face_bc{};

        // Optional, only for CellType::POLYHEDRON
        std::vector<id_t> poly_face_offsets{}; ///< size = n_faces + 1, offsets[0]=0
        std::vector<id_t> poly_face_nodes{};   ///< flattened global node ids of faces
    };

    struct DraftData {
        std::vector<GridNode::Ptr> nodes{};
        std::vector<DraftCell>     cells{};

        /// @brief Clear all draft buffers
        void clear();
    };

    // --------------------------
    // Grid Modification
    // --------------------------

    /// @brief Достроить зеркальное отражение сетки. Все узлы сетки должны
    /// лежать в одной полуплоскости.
    /// @param axis = 0, 1, 2 -- нормали X, Y, Z
    void mirror_(int axis);

    /// @brief Триангуляция Делоне для сетки из четырёхугольников
    void triangulation_quad_delaunay_();

    /// @brief Триангуляция с центральным узлом для сетки из четырёхугольников
    void triangulation_quad_symmetry_();

    /// @brief Преобразовать сетку из четырёхугольников в AMR сетку
    void make_amr_2D_();

    /// @brief Преобразовать сетку из шестигранников в AMR сетку
    void make_amr_3D_();

    /// @brief Преобразовать двумерную AMR сетку в сетку из четырёхугольников
    void reduce_amr_2D_();

    /// @brief Преобразовать трёхмерную AMR сетку в сетку из шестигранников
    void reduce_amr_3D_();

    // ==========================
    // Internal finalize steps
    // ==========================

    /// @brief Ensure editable state; throws if finalized.
    void require_editable_() const;

    /// @brief Ensure finalized state; throws if editable.
    void require_finalized_() const;

    /// @brief Clear draft and all optional data
    void remove_options_();

private:
    State m_state{State::editable};

    std::optional<DraftData> m_draft{std::nullopt};

    // Finalized storage
    int m_dim{0};
    Type m_type{Type::NONE};
    BuildOptions m_options;
    std::vector<Node> m_nodes{};
    std::vector<Cell> m_cells{};

    // Optional caches
    std::optional<GeometryCache> m_geom{std::nullopt};
    std::optional<FacesCache> m_faces{std::nullopt};
    std::optional<EdgesCache> m_edges{std::nullopt};
    std::optional<IncidenceCache> m_inc{std::nullopt};
};

} // namespace zephyr::geom