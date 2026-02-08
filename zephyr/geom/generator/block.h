#pragma once

#include <vector>
#include <memory>
#include <sys/stat.h>

#include <zephyr/mesh/side.h>
#include <zephyr/geom/generator/base_node.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>

namespace zephyr::geom::generator {

using mesh::Side2D;

constexpr int to_axis(Side2D side) {
    return (side == Side2D::B || side == Side2D::T) ? 0 : 1;
}

constexpr std::array<Side2D, 2> sides(int axis) {
    return axis == 0? std::array{Side2D::B, Side2D::T} : std::array{Side2D::L, Side2D::R};
}

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

/// @brief Двумерный массив с возможностью индексации отрицательными
/// числами как в python.
template <typename T>
class Array2D {
public:
    /// @brief Пустой массив
    Array2D() = default;

    /// @brief Создать массив заданного размера
    explicit Array2D(const std::array<int, 2> sizes, const T& value = {})
        : m_sizes(sizes) {
        z_assert(sizes[0] > 0 && sizes[1] > 0, "Bad sizes constructor");
        m_data.resize(m_sizes[0] * m_sizes[1], value);
    }

    /// @brief Пустой массив?
    bool empty() const { return m_data.empty(); }

    /// @brief Изменить размер массива
    void resize(const std::array<int, 2> sizes, const T& value = {}) {
        z_assert(sizes[0] > 0 && sizes[1] > 0, "Bad sizes resize");
        m_sizes = sizes;
        m_data.resize(m_sizes[0] * m_sizes[1], value);
    }

    /// @brief Размер массива по первой оси
    int size1() const { return m_sizes[0]; }

    /// @brief Размер массива по второй оси
    int size2() const { return m_sizes[1]; }

    /// @brief Размер массива по оси axis
    int size(int axis) const { return m_sizes[axis]; }

    /// @brief Получить значение по ссылке
    T& operator()(int i, int j) {
        z_assert(-m_sizes[0] <= i && i < m_sizes[0], "Bad index i");
        z_assert(-m_sizes[1] <= j && j < m_sizes[1], "Bad index j");
        return m_data[reduce(j, m_sizes[1]) * m_sizes[0] + reduce(i, m_sizes[0])];
    }

    /// @brief Получить значение
    const T& operator()(int i, int j) const {
        z_assert(-m_sizes[0] <= i && i < m_sizes[0], "Bad index i");
        z_assert(-m_sizes[1] <= j && j < m_sizes[1], "Bad index j");
        return m_data[reduce(j, m_sizes[1]) * m_sizes[0] + reduce(i, m_sizes[0])];
    }

private:
    /// @brief Привести индекс к корректному диапазону
    static int reduce(int idx, int size) { return idx < 0 ? idx + size : idx; }

    std::array<int, 2> m_sizes{}; ///< Размеры по двум осям
    std::vector<T> m_data{};      ///< Массив с данными
};

/// @brief Представление четырехугольного блока.
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

    /// @brief Число ячеек вдоль граней (v1, v2) и (v3, v4)
    int size1() const { return m_sizes[0]; }

    /// @brief Число ячеек вдоль граней (v1, v3) и (v2, v4)
    int size2() const { return m_sizes[1]; }

    /// @brief Число ячеек по оси блока
    int size(int axis) const { return m_sizes[axis]; }

    /// @brief Число ячеек вдоль стороны
    int size(Side2D side) const { return m_sizes[to_axis(side)]; }

    /// @brief Число ячеек вдоль грани
    int size(BaseNode::Ref v1, BaseNode::Ref v2) const {
        return m_sizes[get_axis(v1, v2)];
    }

    /// @brief Установить число ячеек вдоль оси
    void set_size(int axis, int N);

    /// @brief Установить число ячеек вдоль грани
    void set_size(Side2D side, int N);

    /// @brief Установить число ячеек вдоль грани (v1, v2)
    void set_size(BaseNode::Ref v1, BaseNode::Ref v2, int N);


    /// @brief Установить границу на сторону
    void set_boundary(Side2D side, Curve::Ref curve);

    /// @brief Установить границу на грань (v1, v2)
    void set_boundary(BaseNode::Ref v1, BaseNode::Ref v2, Curve::Ref curve);


    /// @brief Центр блока (среднее базисных вершин)
    Vector3d center() const;

    /// @brief Центр стороны (проецируется на границу, если есть)
    Vector3d center(Side2D side) const;

    /// @brief Получить вершину блока
    BsVertex::Ptr& vertex(int i, int j) { return m_vertices(i, j); }
    
    /// @brief Получить вершину блока
    BsVertex::Ref vertex(int i, int j) const { return m_vertices(i, j); }

    /// @brief Ссылка на двумерный массив вершин
    const Array2D<BsVertex::Ptr>& vertices() const { return m_vertices; }

    /// @brief Стороны блока, которые прилегают к заданной вершине.
    /// Гарантируется перечисление против часовой стрелки (внутри блока)
    std::tuple<Side2D, Side2D> adjacent_sides(const BaseNode* v) const;

    /// @brief Стороны блока, которые прилегают к заданной вершине
    std::tuple<Side2D, Side2D> adjacent_sides(BaseNode::Ref v) const;

    /// @brief Пара прилегающих вершин.
    /// Гарантируется перечисление против часовой стрелки (внутри блока)
    std::tuple<BaseNode::Ptr, BaseNode::Ptr> adjacent_nodes(const BaseNode* v) const;

    /// @brief Внешняя граница ячейки?
    /// @param side Сторона ячейки
    bool is_boundary(Side2D side) const { return m_boundaries[side] != nullptr; }


    double modulus() const { return m_modulus; }

    double* ratio_ptr(int axis) { return &m_ratio[axis]; }

    double* ratio_ptr(Side2D side) { return &m_ratio[to_axis(side)]; };

    double ratio() const { return static_cast<double>(size2()) / static_cast<double>(size1()); }

    /// @brief Угол на границе?
    bool is_boundary(const BaseNode* v) const;

    /// @brief Индекс базисной вершины внутри блока
    int base_node_index(const BaseNode* v) const;

    /// @brief Индекс базисной вершины внутри блока
    int base_node_index(BaseNode::Ref v) const { return base_node_index(v.get()); }

    double calc_modulus() const;

    /// @brief Пересчитывает modulus и ratio1, ratio2
    void update_modulus();


    int index;

    /// @brief Управляющий класс
    friend class BlockStructured;

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


    int get_axis(BaseNode::Ref v1, BaseNode::Ref v2) const;

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
    BsVertex::Ptr& corner_vertex(int v_idx);

    /// @brief Граничная вершина блока
    BsVertex::Ptr& boundary_vertex(Side2D side, int idx);

    /// @brief Вершина со второго ряда от границы
    BsVertex::Ptr& preboundary_vertex(Side2D side, int idx);

private:

    // -------- Основные топологические параметры --------

    /// @brief Базовые вершины (Z-ordering)
    std::array<BaseNode::Ptr, 4> m_base_nodes{}; // nullptr

    /// @brief Границы области (nullptr для внутренних границ)
    std::array<Curve::Ptr, 4> m_boundaries{};

    /// @brief Ссылки на соседние блоки (nullptr для границ области)
    std::array<Block::WPtr, 4> m_adjacent_blocks{};

    /// @brief Поворот соседнего блока [0..3]
    std::array<int, 4> m_rotations{};

    /// @brief Конформный модуль блока (геометрическая характеристика)
    double m_modulus{1.0};

    // Число ячеек вдоль граней (v1, v2) и (v3, v4)
    // Число ячеек вдоль вершин (v1, v3) и (v2, v4)
    std::array<int, 2> m_sizes{0, 0};

    /// @brief Коэффициенты для анизотропного сглаживания Лапласа
    /// modulus * (size2 / size1)
    /// (size1 / size2) / modulus
    std::array<double, 2> m_ratio{1.0, 1.0};

    /// @brief Двумерный массив с вершинами
    Array2D<BsVertex::Ptr> m_vertices{};
};

} // namespace zephyr::mesh::generator
