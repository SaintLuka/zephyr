#pragma once

#include <memory>

#include <zephyr/geom/generator/array2d.h>
#include <zephyr/geom/generator/base_node.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>

namespace zephyr::geom::generator {

class BlockStructured;

/// @brief Пара величин, присвоенная осям блока
template <typename T>
class AxisPair {
public:
    /// @brief Конструктор по умолчанию
    AxisPair() = default;

    /// @brief Конструктор от initializer_list
    constexpr AxisPair(std::initializer_list<T> list) {
        if (list.size() != 2)
            throw std::invalid_argument("Pair requires exactly 2 elements");
        std::copy_n(list.begin(), 2, values.begin());
    }

    /// @brief Конструктор от std::array
    constexpr AxisPair(const std::array<T, 2>& arr) : values(arr) {}

    /// @brief Неявное приведение к std::array
    operator const std::array<T, 2>&() const { return values; }

    /// @brief Доступ по оси
    T& operator[](Axis axis) {
        return values[static_cast<int>(axis)];
    }

    /// @brief Доступ по оси
    const T& operator[](Axis axis) const {
        return values[static_cast<int>(axis)];
    }

    /// @brief Доступ по стороне
    T& operator[](Side side) {
        return values[static_cast<int>(to_axis(side))];
    }

    /// @brief Доступ по стороне
    const T& operator[](Side side) const {
        return values[static_cast<int>(to_axis(side))];
    }

private:
    /// @brief Пара значений
    std::array<T, 2> values{};
};

/// @brief Представление четырехугольного блока.
class Block {
public:
    using Ptr = std::shared_ptr<Block>;
    using Ref = const std::shared_ptr<Block>&;
    using WPtr = std::weak_ptr<Block>;

    /// @brief Конструктор структурированного блока.
    /// @details Вершины сортируются, чтобы получить обход против часовой
    /// стрелки, первая вершина остается на месте.
    explicit Block(const std::array<BaseNode::Ptr, 4>& vertices);

    /// @brief Инициализация структурированного блока.
    /// @details Вершины сортируются, чтобы получить обход против часовой
    /// стрелки, первая вершина остается на месте.
    static Block::Ptr create(const std::array<BaseNode::Ptr, 4>& vertices);

    /// @brief Индекс блока
    int index() const { return m_index; }

    /// @brief Установить индекс блока
    void set_index(int i);

    // ------------------------------ Геометрия -------------------------------

    /// @brief Центр блока (среднее базисных вершин)
    Vector3d center() const;

    /// @brief Центр стороны (проецируется на границу, если есть)
    Vector3d center(Side side) const;

    /// @brief Длина стороны блока
    double length(Side side) const;

    // --------------------------- Базисные вершины ---------------------------

    /// @brief Индекс базисной вершины внутри блока
    int base_node_index(const BaseNode* v) const;

    /// @brief Индекс базисной вершины внутри блока
    int base_node_index(BaseNode::Ref v) const { return base_node_index(v.get()); }

    /// @brief Стороны блока, которые прилегают к заданной вершине.
    /// Гарантируется перечисление против часовой стрелки (внутри блока)
    std::tuple<Side, Side> adjacent_sides(const BaseNode* v) const;

    /// @brief Стороны блока, которые прилегают к заданной вершине.
    std::tuple<Side, Side> adjacent_sides(BaseNode::Ref v) const;

    /// @brief Пара прилегающих вершин. Гарантируется перечисление вершин против
    /// часовой стрелки (внутри блока).
    std::tuple<BaseNode::Ptr, BaseNode::Ptr> adjacent_nodes(const BaseNode* v) const;

    /// @brief Базовые вершины
    const auto& base_nodes() const { return m_base_nodes; }

    /// @brief Базовая вершина
    /// @param v_idx Индекс вершины
    BaseNode::Ptr base_node(int v_idx) const { return m_base_nodes[v_idx]; }

    /// @brief Базисная вершина по стороне (против часовой стрелки)
    BaseNode::Ptr base_node(Side side, int idx) const;


    // ---------------------------- Границы блока -----------------------------

    /// @brief Кривая границы области (допускается nullptr)
    Curve::Ref boundary(Side side) const;

    /// @brief Сторона принадлежит границе области?
    bool is_boundary(Side side) const;

    /// @brief Угол принадлежит границе области?
    bool is_boundary(const BaseNode* v) const;

    /// @brief Установить границу на сторону
    void set_boundary(Side side, Curve::Ref curve);

    /// @brief Установить границу на сторону (v1, v2)
    void set_boundary(BaseNode::Ref v1, BaseNode::Ref v2, Curve::Ref curve);


    // --------------------- Выбор сторон и смежные блоки ---------------------

    /// @brief Индекс грани
    /// @param v1, v2 Вершины грани
    Side get_side(BaseNode::Ref v1, BaseNode::Ref v2) const;

    /// @brief Ось выбранной грани
    /// @param v1, v2 Вершины грани
    Axis get_axis(BaseNode::Ref v1, BaseNode::Ref v2) const;

    /// @brief Смежный блок, может вернуть nullptr (для expired).
    /// @param side Сторона ячейки
    Block::Ptr adjacent_block(Side side) const;

    /// @brief Сторона этой же грани у соседнего блока
    Side twin_face(Side side) const;


    // ----------------------- Функции "верхнего уровня" ------------------------

    /// @brief Связать два блока
    static void link(Block::Ref B1, Block::Ref B2);

    /// @brief Сгенерировать узлы сетки с нуля
    Array2D<BsVertex::Ptr> create_vertices(AxisPair<int> sizes) const;

    /// @brief Сгенерировать узлы сетки на основе сохраненного отображения
    Array2D<BsVertex::Ptr> create_vertices_again(AxisPair<int> sizes) const;


    // -------------------------- Конформные приколы --------------------------

    /// @brief Конформный модуль криволинейного четырехугольника
    double modulus() const { return m_modulus; }

    /// @brief Установить конформный модуль (пересчитать lambda)
    void set_modulus(double K);

    /// @brief Оценить и установить конформный модуль четырехугольника
    double estimate_modulus() const;

    /// @brief Пересчитывает модуль и коэффициенты сглаживания
    void update_modulus(const Array2D<BsVertex::Ptr>& vertices);

    /// @brief Получить конформное отображение блока
    const Array2D<Vector3d>& mapping() const { return m_mapping; }

    /// @brief Сохранить текущее отображение
    void set_mapping(const Array2D<BsVertex::Ptr>& vertices);

private:
    // Основные топологические параметры, задаются пользователем

    /// @brief Индекс блока в общем списке блоков
    int m_index{-1};

    /// @brief Базовые вершины (Z-ordering)
    std::array<BaseNode::Ptr, 4> m_base_nodes{};

    /// @brief Границы области (nullptr для внутренних границ)
    std::array<Curve::Ptr, 4> m_boundaries{};

    // Связи блоков, находятся при вызове link()

    /// @brief Ссылки на соседние блоки (nullptr для границ области)
    std::array<Block::WPtr, 4> m_adjacent_blocks{};

    /// @brief Поворот соседнего блока [0..3]
    std::array<int, 4> m_rotations{};

    /// @brief Конформный модуль блока (геометрическая характеристика)
    double m_modulus{NAN};

    /// @brief Конформное отображение блока в виде таблицы вершин
    Array2D<Vector3d> m_mapping{};
};

} // namespace zephyr::mesh::generator
