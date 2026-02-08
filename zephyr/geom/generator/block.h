#pragma once

#include <vector>
#include <memory>

#include <zephyr/mesh/side.h>
#include <zephyr/geom/generator/base_node.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>

namespace zephyr::geom::generator {

using mesh::Side2D;

constexpr int to_face_idx(Side2D side) {
    switch (side) {
        case Side2D::BOTTOM: return 0;
        case Side2D::RIGHT: return 1;
        case Side2D::TOP: return 2;
        case Side2D::LEFT: return 3;
        default: throw std::invalid_argument("Invalid side");
    }
}

constexpr Side2D idx_to_side(int idx) {
    switch (idx) {
        case 0: return Side2D::BOTTOM;
        case 1: return Side2D::RIGHT;
        case 2: return Side2D::TOP;
        case 3: return Side2D::LEFT;
        default: throw std::invalid_argument("Invalid side");
    }
}

inline int mod(int idx, int size) {
    return idx < 0 ? idx + size : idx;
}

class RotatedTable2D;

template <typename T>
class Table2D {
public:
    Table2D() = default;

    Table2D(int size1, int size2, const T& value = {})
        : m_size1(size1), m_size2(size2) {
        m_data.resize(m_size1 * m_size2, value);
    }

    void resize(int size1, int size2, const T& value = {}) {
        m_size1 = size1;
        m_size2 = size2;
        m_data.resize(m_size1 * m_size2, value);
    }


    int size1() const { return m_size1; }
    int size2() const { return m_size2; }

    bool empty() const { return m_data.empty(); }

    T& get(int i, int j) { return m_data[mod(j, m_size2) * m_size1 + mod(i, m_size1)]; }

    const T& get(int i, int j) const { return m_data[mod(j, m_size2) * m_size1 + mod(i, m_size1)]; }

    T& operator()(int i, int j) { return get(i, j); }

    const T& operator()(int i, int j) const { return get(i, j); }


    void print() const {
        for (int j = m_size2 - 1; j >= 0; --j) {
            for (int i = 0; i < m_size1; ++i) {
                std::cout << get(i, j) << " ";
            }
            std::cout << "\n";
        }
    }

    // Proxy-object (view)
    class Rotated {
    public:
        Table2D& ref;
        int r;

        Rotated(Table2D& ref, int r) : ref(ref), r(r) {}

        int size1() const { return r % 2 ? ref.size2() : ref.size1(); }
        int size2() const { return r % 2 ? ref.size1() : ref.size2(); }

        T& get(int i, int j) {
            switch (mod(r, 4)) {
                case 1: return ref.get(j, ref.size2() - 1 - i);
                case 2: return ref.get(ref.size1() - 1 - i, ref.size2() - 1 - j);
                case 3: return ref.get(ref.size1() - 1 - j, i);
                default: return ref.get(i, j);
            }
        }

        const T& get(int i, int j) const {
            switch (mod(r, 4)) {
                case 1: return ref.get(j, ref.size2() - 1 - i);
                case 2: return ref.get(ref.size1() - 1 - i, ref.size2() - 1 - j);
                case 3: return ref.get(ref.size1() - 1 - j, i);
                default: return ref.get(i, j);
            }
        }

        T& operator()(int i, int j) { return get(i, j); }

        const T& operator()(int i, int j) const { return get(i, j); }

        void print() const {
            for (int j = size2() - 1; j >= 0; --j) {
                for (int i = 0; i < size1(); ++i) {
                    std::cout << get(i, j) << " ";
                }
                std::cout << "\n";
            }
        }
    };

    Rotated R(int r) { return Rotated{*this, r}; }

private:
    int m_size1, m_size2;
    std::vector<T> m_data;
};

/// @brief Представление четырехугольного блока.
/// Основа блочно-структурированной сетки, большинство функций класса открыты
/// только для управляющего класса BlockStructured.
class Block : public std::enable_shared_from_this<Block> {
public:
    using Ptr = std::shared_ptr<Block>;
    using Ref = const std::shared_ptr<Block>&;
    using WPtr = std::weak_ptr<Block>;

    /// @brief Пустой блок
    Block() = default;

    /// @brief Инициализация структурированного блока.
    /// @details Вершины сортируются, чтобы получить обход против часовой
    /// стрелки, первая вершина остается на месте.
    Block &operator=(const std::array<BaseNode::Ptr, 4>& vertices);

    /// @brief Число ячеек вдоль грани
    /// @param v1, v2 Вершины грани (v1, v2)
    int size(BaseNode::Ref v1, BaseNode::Ref v2) const;

    /// @brief Установить число ячеек вдоль грани (v1, v2)
    void set_size(BaseNode::Ref v1, BaseNode::Ref v2, int N);

    /// @brief Установить границу на грань (v1, v2)
    void set_boundary(BaseNode::Ref v1, BaseNode::Ref v2, Curve::Ref curve);

    /// @brief Центр блока
    Vector3d center() const;

    /// @brief Центр стороны (проецируется на границу, если есть)
    Vector3d center(Side2D side) const;
    
    /// @brief Получить вершину
    BsVertex::Ptr vertex(int i, int j) const;

    /// @brief Пара прилегающих вершин.
    /// Гарантируется перечисление против часовой стрелки (внутри блока)
    std::tuple<BaseNode::Ptr, BaseNode::Ptr> adjacent_nodes(const BaseNode* v) const;

    /// @brief Стороны блока, которые прилегают к заданной вершине.
    /// Гарантируется перечисление против часовой стрелки (внутри блока)
    std::tuple<Side2D, Side2D> adjacent_sides(const BaseNode* v) const;

    /// @brief Стороны блока, которые прилегают к заданной вершине
    std::tuple<Side2D, Side2D> adjacent_sides(BaseNode::Ref v) const;

    /// @brief Внешняя граница ячейки?
    /// @param side Сторона ячейки
    bool is_boundary(Side2D side) const { return m_boundaries[side] != nullptr; }


    double modulus() const { return m_modulus; }

    double ratio() const { return static_cast<double>(size2()) / static_cast<double>(size1()); }

    /// @brief Угол на границе?
    bool is_boundary(const BaseNode* v) const;

    /// @brief Индекс базисной вершины внутри блока
    int node_index(const BaseNode* v) const;

    /// @brief Индекс базисной вершины внутри блока
    int node_index(BaseNode::Ref v) const { return node_index(v.get()); }

    double calc_modulus() const;


    int index;

    /// @brief Управляющий класс
    friend class BlockStructured;

    /// @brief Разбиение вдоль граней (v1, v2) и (v3, v4)
    int size1() const { return m_size1; }

    /// @brief Разбиение вдоль граней (v1, v4) и (v2, v3)
    int size2() const { return m_size2; }

    /// @brief Число ячеек вдоль стороны
    /// @param side Сторона ячейки
    int get_size(Side2D side) const;

    /// @brief Установить число ячеек вдоль грани
    void set_size(Side2D side, int N);

    /// @brief Смежный блок, может вернуть nullptr (для expired).
    /// @param side Сторона ячейки
    Block::Ptr adjacent_block(Side2D side) const;

    /// @brief Граница ячейки
    /// @param side Сторона ячейки
    Curve::Ref boundary(Side2D side) const { return m_boundaries[side]; }

    /// @brief Базовые вершины
    const auto& base_nodes() const { return m_base_nodes; }

    /// @brief Базовая вершина
    /// @param v_idx Индекс вершины
    BaseNode::Ptr base_node(int v_idx) const { return m_base_nodes[v_idx]; }

    /// @brief Индекс грани
    /// @param v1, v2 Вершины грани
    Side2D get_side(BaseNode::Ref v1, BaseNode::Ref v2) const;

    /// @brief Сторона такой же грани у соседа
    Side2D neib_face(Side2D side) const;

    /// @brief Связать два блока
    static void link(Block::Ref B1, Block::Ref B2);

    /// @brief Сгенерировать узлы сетки
    void create_vertices();

    /// @brief Связать узлы сетки
    void link_vertices();

    /// @brief Сглаживание вершин в блоке
    /// @return Максимальный относительный сдвиг вершин
    double smooth();

    /// @brief Обновить положение вершин
    void update();

    /// @brief Угловая вершина блока
    BsVertex::Ptr &corner_vertex(int v_idx);

    /// @brief Граничная вершина блока
    BsVertex::Ptr &boundary_vertex(Side2D side, int idx);

    /// @brief Вершина со второго ряда от границы
    BsVertex::Ptr &preboundary_vertex(Side2D side, int idx);

private:

    // -------- Основные топологические параметры --------

    /// @brief Базовые вершины обходятся против часовой стрелки
    std::array<BaseNode::Ptr, 4> m_base_nodes{}; // nullptr

    /// @brief Границы области (nullptr для внутренних границ)
    std::array<Curve::Ptr, 4> m_boundaries{}; // nullptr

    /// @brief Ссылки на соседние блоки (nullptr для границ области)
    std::array<Block::WPtr, 4> m_adjacent_blocks{}; // nullptr

    /// @brief Поворот соседнего блока [0..3]
    std::array<int, 4> m_rotations{-1, -1, -1, -1};

    // ???
    double m_modulus{1.0};


    // -------- Дополнительные параметры --------

    int m_size1{0};  ///< Число ячеек вдоль граней (v1, v2) и (v3, v4)
    int m_size2{0};  ///< Число ячеек вдоль вершин (v1, v4) и (v2, v3)

    /// @brief Двумерный массив с вершинами
    Table2D<BsVertex::Ptr> m_vertices{};
};

} // namespace zephyr::mesh::generator
