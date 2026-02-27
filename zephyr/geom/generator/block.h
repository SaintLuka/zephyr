#pragma once

#include <memory>

#include <zephyr/geom/generator/array2d.h>
#include <zephyr/geom/generator/base_node.h>
#include <zephyr/geom/generator/curve/curve.h>

namespace zephyr::geom::generator {

class BsVertex;

/// @brief Таблица указателей на внутренние узлы
using Table2D = Array2D<std::shared_ptr<BsVertex>>;

/// @brief Четырёхугольный блок на базисных вершинах
class Block {
    using BsVertex_Ptr = std::shared_ptr<BsVertex>;
public:
    using Ptr = std::shared_ptr<Block>;
    using Ref = const std::shared_ptr<Block>&;
    using WPtr = std::weak_ptr<Block>;

    /// @brief Конструктор структурированного блока.
    /// @param nodes Базисные вершины блока в произвольном порядке. Внутри
    /// конструктора упорядочиваются, первая вершина остается на месте.
    explicit Block(const std::array<BaseNode::Ptr, 4>& nodes);

    /// @brief Конструктор структурированного блока.
    /// @param nodes Базисные вершины блока в произвольном порядке. Внутри
    /// конструктора упорядочиваются, первая вершина остается на месте.
    static Block::Ptr create(const std::array<BaseNode::Ptr, 4>& nodes);

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

    /// @brief Базовые вершины
    const auto& base_nodes() const { return m_base_nodes; }

    /// @brief Базовая вершина
    /// @param v_idx Индекс вершины
    BaseNode::Ptr base_node(int v_idx) const { return m_base_nodes[v_idx]; }

    /// @brief Базисная вершина по стороне (индексация против часовой стрелки)
    BaseNode::Ptr base_node(Side side, int idx) const;

    /// @brief Индекс базисной вершины внутри блока
    int base_node_index(const BaseNode* v) const;

    /// @brief Индекс базисной вершины внутри блока
    int base_node_index(BaseNode::Ref v) const { return base_node_index(v.get()); }

    /// @brief Стороны блока, которые прилегают к заданной вершине.
    /// Гарантируется перечисление против часовой стрелки (внутри блока).
    std::tuple<Side, Side> incident_sides(const BaseNode* v) const;

    /// @brief Стороны блока, которые прилегают к заданной вершине.
    /// Гарантируется перечисление против часовой стрелки (внутри блока).
    std::tuple<Side, Side> incident_sides(BaseNode::Ref v) const;

    /// @brief Пара прилегающих вершин. Гарантируется перечисление вершин против
    /// часовой стрелки (внутри блока).
    std::tuple<BaseNode::Ptr, BaseNode::Ptr> adjacent_nodes(const BaseNode* v) const;

    /// @brief Пара прилегающих вершин. Гарантируется перечисление вершин против
    /// часовой стрелки (внутри блока).
    std::tuple<BaseNode::Ptr, BaseNode::Ptr> adjacent_nodes(BaseNode::Ref v) const {
        return adjacent_nodes(v.get());
    }

    // ---------------------------- Границы блока -----------------------------

    /// @brief Кривая границы области. Возвращает nullptr при отсутствии
    /// границы, поэтому может использоваться в условиях.
    Curve::Ref boundary(Side side) const;

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

    // ---------------------- Функции "верхнего уровня" -----------------------

    /// @brief Связать два блока
    static void link(Block::Ref B1, Block::Ref B2);

    /// @brief Сгенерировать таблицу узлов сетки, если в блоке определено
    /// конформное отображение, то оно будет использовано для генерации.
    Table2D create_vertices(AxisPair<int> sizes) const;

    // -------------------------- Конформные приколы --------------------------

    /// @brief Конформный модуль криволинейного четырехугольника
    double modulus() const { return m_modulus; }

    /// @brief Оценить и установить конформный модуль четырехугольника
    void estimate_modulus();

    /// @brief Пересчитывает модуль и коэффициенты сглаживания
    void update_modulus(const Table2D& vertices);

    /// @brief Сохранить текущее отображение
    void set_mapping(const Table2D& vertices);

    /// @brief Получить конформное отображение блока
    const Array2D<Vector3d>& mapping() const { return m_mapping; }

private:
    /// @brief Установить конформный модуль
    void set_modulus(double K);

    /// @brief Сгенерировать узлы сетки с нуля
    Table2D create_vertices_init(AxisPair<int> sizes) const;

    /// @brief Сгенерировать узлы сетки на основе сохраненного отображения
    Table2D create_vertices_again(AxisPair<int> sizes) const;

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

/// @brief Пара смежных блоков
struct BlockPair {
    /// @brief Пара смежных блоков
    Block::WPtr b1{}, b2{};

    /// @brief Стороны у смежных блоков
    Side side1{}, side2{};

    /// @brief Смежные блоки не заданы
    bool empty() const { return b1.expired(); }

    /// @brief У ребра только один блок (граница)
    bool boundary() const { return b2.expired(); }

    /// @brief У ребра два блока (внутреннее ребро)
    bool inner() const { return b2.expired(); }

    /// @brief Добавить блок (не кидает ошибок)
    void add(Block::Ref block, Side side);
};

} // namespace zephyr::mesh::generator
