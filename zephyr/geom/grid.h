#pragma once

#include <memory>
#include <vector>
#include <optional>
#include <boost/container/static_vector.hpp>

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

/// @brief Вектор неопределенных значений
inline Vector3d nanvec() { return {NAN, NAN, NAN}; }

/// @brief Тип локального индекса в ячейке
using node_id_t = std::uint8_t;

/// @brief Максимальное число узлов в ячейке (до 255)
constexpr std::uint8_t max_nodes_per_cell = std::numeric_limits<std::uint8_t>::max();

/// @brief Максимальное число вершин на грани ячейки
constexpr std::uint8_t max_nodes_per_face = 16;

/// @brief Локальные индексы вершин внутри ячейки (можно использовать std::vector)
using node_ids_t = boost::container::static_vector<node_id_t, max_nodes_per_face>;

// ============================================================================
//                           Сеточные примитивы
// ============================================================================

/// @brief Узел сетки (положение + индекс). Индекс проставляет Grid,
/// при этом выполняется дедубликация по указателям.
class Node {
public:
    using Ptr = std::shared_ptr<Node>;
    using Ref = const std::shared_ptr<Node>&;

    Vector3d pos{nanvec()};           ///< Положение
    Boundary bc{Boundary::UNDEFINED}; ///< Граничное условие


    /// @brief Создать узел из Eigen вектора
    explicit Node(const Vector3d& v);

    /// @brief Создать узел из Eigen вектора
    static Ptr create(const Vector3d& v);

    /// @brief Создать узел (x, y) на плоскости
    static Ptr create(double x, double y);

    /// @brief Индекс узла
    id_t id() const noexcept { return m_id; }

    /// @brief Конечные значения? Нет NAN и infinity.
    bool is_finite() const;

    /// @brief Точное сравнение координат и индексов
    bool operator==(const Node& other) const;

    /// @brief Точное сравнение координат и индексов
    bool operator!=(const Node& other) const;

private:
    friend class Grid;
    id_t m_id{invalid_id};  ///< Индекс (для Grid builder)
};

/// @brief Грань ячейки
class Face {
public:
    /// @brief Пустая грань без узлов и гран условий
    Face() = default;

    // ------------------------------- Индексы --------------------------------

    /// @brief Локальные индексы проставлены?
    bool has_nodes() const { return !m_nodes.empty(); }

    /// @brief Число узлов грани
    int n_nodes() const { return static_cast<int>(m_nodes.size()); }

    /// @brief Получить локальный индекс i-ого узла
    int node_idx(int i) const { return m_nodes[i]; }

    /// @brief Локальные индексы узлов грани
    const node_ids_t& nodes() const { return m_nodes; }

    /// @brief Установить локальные индексы узлов
    void set_nodes(std::span<const int> node_ids);

    // ------------------------------ Связность -------------------------------

    /// @brief Граничное условие на грани
    Boundary bc() const { return m_bc; }

    /// @brief Установить граничное условие
    void set_bc(Boundary bc);

    /// @brief Индекс соседней ячейки
    id_t neib() const { return m_neib; }

    /// @brief Установить индекс соседней ячейки (neib_id)
    /// и индекс смежной грани внутри этой ячейки (face_id).
    void set_neib(id_t neib_id, int face_id);

    // ------------------------------ Геометрия -------------------------------

    /// @brief Площадь грани
    double area() const { return m_area_n.norm(); }

    /// @brief Внешняя нормаль (относительно ячейки-владельца)
    Vector3d normal() const { return m_area_n.normalized(); }

    /// @brief Барицентр грани
    Vector3d center() const { return m_center; }

    /// @brief Установить площадь и внешнюю нормаль
    void set_area_n(const Vector3d& S) { m_area_n = S; }

    /// @brief Установить центр ячейки
    void set_center(const Vector3d& C) { m_center = C; }

private:
    node_ids_t m_nodes{};           ///< Локальные индексы (внутри ячейки)
    Vector3d   m_area_n{nanvec()};  ///< Площадь на внешнюю нормаль
    Vector3d   m_center{nanvec()};  ///< Барицентр грани

    Boundary m_bc{Boundary::UNDEFINED};  ///< Граничное условие
    id_t m_neib{invalid_id};             ///< Индекс соседней ячейки (-1 для boundary)
    int  m_twin{-1};                     ///< Индекс грани у соседней ячейки
};

/// @brief Ячейка сетки.
class Cell {
public:
    // ---------------------------- Инициализация -----------------------------

    /// @brief Конструктор базовой ячейки (только необходимое)
    Cell(CellType type, std::vector<id_t>&& node_ids);

    /// @brief Конструктор базовой ячейки (только необходимое)
    Cell(CellType type, std::initializer_list<id_t> node_ids);

    /// @brief Установить граничные условия, если граней нет, то они будут созданы
    void set_face_bc(const std::vector<Boundary>& face_bc);

    /// @brief Инициализировать локальные индексы вершин граней (не работает для POLYHEDRON)
    void init_faces();

    /// @brief Установить локальные индексы вершин граней (только для POLYHEDRON)
    void set_faces(const std::vector<std::vector<int>>& face_ids);


    /// @brief Тип ячейки
    CellType type() const { return m_type; }

    /// @brief Число вершин ячейки
    int n_nodes() const { return static_cast<int>(m_nodes.size()); }

    /// @brief Получить индекс вершины ячейки
    id_t node_idx(int i) const { return m_nodes[i]; }

    /// @brief Массив индексов вершин ячейки
    const std::vector<id_t>& nodes() const { return m_nodes; }

    /// @brief Массив граничных условий константного размера
    template <int N>
    std::array<Boundary, N> faces_bc() const {
        std::array<Boundary, N> bc;
        bc.fill(Boundary::UNDEFINED);
        for (int i = 0; i < std::min(N, n_faces()); ++i) {
            bc[i] = m_faces[i].bc();
        }
        return bc;
    }



    bool has_faces() const { return !m_faces.empty(); }

    int n_faces() const { return static_cast<int>(m_faces.size()); }

    const Face& get_face(int iface) const { return m_faces[iface]; }

    void set_bc(int iface, Boundary bc) { m_faces[iface].set_bc(bc); }

    void set_neib(int iface, id_t neib_id, int face_id);

    const std::vector<Face>& faces() const { return m_faces; }



    // ----------------------------- Модификации ------------------------------

    /// @brief Заменить индексы узлов (количество должно быть то же)
    void replace_nodes(std::vector<id_t>&& new_nodes);

    /// @brief Зеркально отразить ячейку
    void mirror();

    // ------------------------------ Геометрия -------------------------------

    /// @brief Объем ячейки
    double volume() const { return m_volume; }

    /// @brief Барицентр ячейки
    Vector3d centroid() const { return m_center; }

    /// @brief Среднее вершин грани
    Vector3d face_center(const std::vector<Node>& grid_nodes, int iface) const;

    /// @brief Среднее вершин ячейки
    Vector3d center(const std::vector<Node>& grid_nodes) const;

    /// @brief Вычислить геометрию ячейки и граней
    void calc_geom(const std::vector<Node>& grid_nodes);

private:
    CellType m_type{};            ///< Тип ячейки
    std::vector<id_t> m_nodes{};  ///< Индексы вершин
    std::vector<Face> m_faces{};  ///< Грани (optional)

    Vector3d m_center{nanvec()};  ///< Барицентр ячейки
    double   m_volume{NAN};       ///< Объем ячейки
};

// FaceKey for hash-table (order-insensitive: sorted node ids)
struct FaceKey {
    std::vector<id_t> ids;

    explicit FaceKey(const Cell& cell, const Face& face);

    bool operator==(const FaceKey& o) const noexcept;
};

// EdgeKey for hash-table (order-insensitive: sorted node ids)
struct EdgeKey {
    id_t nid1, nid2;

    explicit EdgeKey(id_t i, id_t j);

    bool operator==(const EdgeKey& o) const noexcept;
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
        Editable, Finalized
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
        bool build_faces{false};       ///< Build faces per cell
        bool build_node_cells{false};  ///< Build IncidenceCache::node_cells

        void validate();
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
    bool is_editable() const noexcept { return m_state == State::Editable; }

    /// @brief True if grid is finalized.
    bool is_finalized() const noexcept { return m_state == State::Finalized; }

    // --------------------------
    // Editable API
    // --------------------------

    /// @brief Reset grid to empty editable state (invalidates all data).
    void clear();

    /// @brief Reserve array for nodes
    void reserve_nodes(id_t n_nodes);

    /// @brief Reserve array for cells
    void reserve_cells(id_t n_cells);

    /// @brief Add a node via shared pointer; returns node index.
    id_t add_node(Node::Ref node);

    /// @brief Add a node via shared pointer; returns node indices.
    std::vector<id_t> add_nodes(const std::vector<Node::Ptr>& nodes);

    /// @brief Add a cell using a vector of
    id_t add_cell(CellType type, const std::vector<Node::Ptr>& nodes,
                  const std::vector<Boundary>& faces_bc = {});

    /// @brief Add a cell of special kind (CellType::POLYHEDRON)
    id_t add_polyhedron(const std::vector<Node::Ptr>& nodes,
                        const std::vector<std::vector<Node::Ptr>>& faces,
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
    void triangulation(int mode = 2);

    /// @brief Делает сетку из кубиков сеткой из пирамидок
    void pyramidize();

    /// @brief Транслировать двумерную сетку вдоль направления
    /// @param p Полный вектор сдвига
    /// @param N Число ячеек вдоль ребра
    void extrude(const Vector3d& p, int N, Boundary side1, Boundary side2);

    /// @brief Convert Type::QUAD grid to Type::AMR grid (add additional unique nodes)
    void make_amr();

    // --------------------------
    // Finalize API
    // --------------------------

    /// @brief Build requested extended connectivity/geometry and switch to finalized state.
    ///        Releases draft storage unless keep_input_nodes == true.
    void finalize(const BuildOptions& options);

    // --------------------------
    // Read-only access (finalized)
    // --------------------------

    /// @brief Grid dimension (2 or 3)
    int dimension() const noexcept { return m_dim; }

    /// @brief Adaptive mesh (AMR cells is used)
    bool adaptive() const noexcept { return m_type == Type::AMR; }

    /// @brief 3D polyhedral mesh?
    bool polyhedral() const { return m_type == Type::POLYHEDRON; }

    /// @brief Number of unique nodes (finalized).
    std::size_t n_nodes() const noexcept { return m_nodes.size(); }

    /// @brief Number of cells (finalized).
    std::size_t n_cells() const noexcept { return m_cells.size(); }

    /// @brief Access finalized nodes.
    const std::vector<Node>& nodes() const noexcept { return m_nodes; }

    /// @brief Access finalized cells.
    const std::vector<Cell>& cells() const noexcept { return m_cells; }

    /// @brief Faces is computed?
    bool has_faces() const noexcept;

    /// @brief Number of faces if computed (else 0)
    std::size_t total_faces_per_cell() const noexcept;

    /// @brief Number of nodes for each cell (with duplicates).
    std::size_t total_nodes_per_cell() const noexcept;

private:

    // ==========================
    // Grid Modification
    // ==========================

    /// @brief Создать узлы в центрах ячеек и добавить их на сетку
    std::vector<Node::Ptr> create_central_nodes_();

    /// @brief Создать узлы в центрах граней и добавить их на сетку
    std::unordered_map<FaceKey, Node::Ptr> create_face_nodes_();

    /// @brief Создать узлы в центрах рёбер и добавить их на сетку
    std::unordered_map<EdgeKey, Node::Ptr> create_edge_nodes_();

    /// @brief Достроить зеркальное отражение сетки. Все узлы сетки должны
    /// лежать в одной полуплоскости.
    /// @param axis = 0, 1, 2 -- нормали X, Y, Z
    void mirror_(int axis);

    /// @brief Триангуляция Делоне для сетки из четырёхугольников
    void triangulation_quad_delaunay_();

    /// @brief Триангуляция с центральным узлом для сетки из четырёхугольников
    void triangulation_quad_symmetry_();

    /// @brief Триангуляция с центральным узлом для сетки из шестигранников
    void triangulation_hex_symmetry_();

    /// @brief Преобразовать сетку из четырёхугольников в AMR сетку
    void make_amr_2D_();

    /// @brief Преобразовать сетку из шестигранников в AMR сетку
    void make_amr_3D_();

    // ==========================
    // Internal finalize steps
    // ==========================

    /// @brief Ensure editable state; throws if finalized.
    void require_editable_() const;

    /// @brief Ensure finalized state; throws if editable.
    void require_finalized_() const;

    void initialize_geom_(const BuildOptions& options);

    void initialize_faces_(const BuildOptions& options);

private:
    State m_state{State::Editable};

    int m_dim{0};
    Type m_type{Type::NONE};
    std::vector<Node> m_nodes{};
    std::vector<Cell> m_cells{};
};

} // namespace zephyr::geom

/// @brief Хэш-функция для FaceKey
template <>
struct std::hash<zephyr::geom::FaceKey> {
    size_t operator()(const zephyr::geom::FaceKey& k) const noexcept {
        size_t h = 1469598103934665603ull;
        for (size_t v: k.ids) {
            h ^= v + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
        }
        return h;
    }
};

/// @brief Хэш-функция для EdgeKey
template <>
struct std::hash<zephyr::geom::EdgeKey> {
    size_t operator()(const zephyr::geom::EdgeKey& k) const noexcept {
        size_t h = 1469598103934665603ull;
        for (size_t v: {k.nid1, k.nid2}) {
            h ^= v + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
        }
        return h;
    }
};