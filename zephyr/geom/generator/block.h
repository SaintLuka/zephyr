#pragma once

#include <vector>
#include <memory>

#include <zephyr/geom/generator/base_node.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>
#include <zephyr/utils/array2d.h>

namespace zephyr::geom::generator {

using utils::Array2D;

/// @brief Две оси четырёхугольного блока
enum class Axis : int { X = 0, Y = 1 };

/// @brief Стороны четырёхугольного блока
enum class Side : int {
    L = 0, LEFT   = 0,
    R = 1, RIGHT  = 1,
    B = 2, BOTTOM = 2,
    T = 3, TOP    = 3,
};

/// @brief Список всех сторон четырёхугольника
constexpr std::array sides_2D = {Side::L, Side::R, Side::B, Side::T};

class BlockStructured;

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

    /// @brief Получить базисные вершины некоторой стороны (порядка нет)
    std::tuple<BaseNode::Ptr, BaseNode::Ptr> base_nodes(Side side) const;


    // ---------------------------- Размеры блока -----------------------------

    /// @brief Число ячеек вдоль сторон (v1, v2) и (v3, v4)
    int size1() const { return m_sizes[0]; }

    /// @brief Число ячеек вдоль сторон (v1, v3) и (v2, v4)
    int size2() const { return m_sizes[1]; }

    /// @brief Число ячеек вдоль оси блока
    int size(Axis axis) const;

    /// @brief Число ячеек вдоль стороны блока
    int size(Side side) const;

    /// @brief Число ячеек вдоль стороны (v1, v2)
    int size(BaseNode::Ref v1, BaseNode::Ref v2) const;

    /// @brief Установить число ячеек вдоль оси блока
    void set_size(Axis axis, int N);

    /// @brief Установить число ячеек вдоль стороны блока
    void set_size(Side side, int N);

    /// @brief Установить число ячеек вдоль стороны (v1, v2)
    void set_size(BaseNode::Ref v1, BaseNode::Ref v2, int N);

    /// @brief Строка информации о разбиении блока на ячейки
    std::string sizes_info() const;

    /// @brief Относительный размер вдоль оси блока
    double rel_size(Axis axis) const;

    /// @brief Относительный размер вдоль стороны (v1, v2)
    double rel_size(BaseNode::Ref v1, BaseNode::Ref v2) const;

    /// @brief Сбросить относительные размеры
    void reset_rel_sizes();

    /// @brief Выставить начальные относительные размеры
    void init_rel_sizes();

    /// @brief Обновить размеры (из конформного модуля)
    /// @return true, если оба размера теперь заданы
    bool update_rel_sizes();

    /// @brief Установить относительный размер вдоль оси блока
    void set_rel_size(Axis axis, double N);

    /// @brief Установить относительный размер вдоль стороны блока
    void set_rel_size(Side side, double N);

    /// @brief Установить относительный размер вдоль стороны (v1, v2)
    void set_rel_size(BaseNode::Ref v1, BaseNode::Ref v2, double N);

    /// @brief Установить относительные размеры блока
    void set_rel_sizes(double Nx, double Ny);


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


    // ----------------------- Внутренние вершины блока -----------------------

    /// @brief Ссылка на двумерный массив вершин
    const Array2D<BsVertex::Ptr>& vertices() const { return m_vertices; }

    /// @brief Получить вершину блока
    BsVertex::Ptr& vertex(int i, int j) { return m_vertices(i, j); }
    
    /// @brief Получить вершину блока
    BsVertex::Ref vertex(int i, int j) const { return m_vertices(i, j); }

    /// @brief Угловая вершина блока
    BsVertex::Ptr& corner_vertex(int v_idx);

    /// @brief Граничная вершина блока
    BsVertex::Ptr& boundary_vertex(Side side, int idx);

    /// @brief Вершина со второго ряда от границы
    BsVertex::Ptr& near_boundary_vertex(Side side, int idx);


    // -------------------------- Конформные приколы --------------------------

    /// @brief Конформный модуль криволинейного четырехугольника
    double modulus() const;

    /// @brief Установить конформный модуль (пересчитать lambda)
    void set_modulus(double K);

    /// @brief Оценить и установить конформный модуль четырехугольника
    void estimate_modulus();

    /// @brief Указатель на коэффициент анизотропного сглаживания
    double* lambda_ptr(Axis axis);

    /// @brief Указатель на коэффициент анизотропного сглаживания
    double* lambda_ptr(Side side);

    /// @brief Обновить коэффициенты сглаживания (по модулю и размерам)
    void update_lambda();

    /// @brief Строка информации о параметрах
    std::string conformal_info() const;


    // ----------------------- Функции "верхнего уровня" ------------------------

    /// @brief Связать два блока
    static void link(Block::Ref B1, Block::Ref B2);

    /// @brief Очистить массив внутренних узлов
    void clear_vertices();

    /// @brief Проверить размеры смежных блоков перед созданием вершин
    void check_consistency() const;

    /// @brief Сгенерировать узлы сетки с нуля
    void create_vertices_init();

    /// @brief Сгенерировать узлы сетки на основе предыдущих
    void create_vertices_again();

    /// @brief Склеить узлы на границах блоков
    void merge_vertices();

    /// @brief Связать узлы сетки
    void link_vertices();

    /// @brief Сглаживание вершин в блоке
    void smooth();

    /// @brief Обновить положение вершин
    /// @return Максимальный относительный сдвиг вершин
    double move_vertices();

    /// @brief Пересчитывает модуль и коэффициенты сглаживания
    void update_modulus();

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
    std::array<Block::WPtr, 4> m_adjacent_blocks{}; // ~nullptr

    /// @brief Поворот соседнего блока [0..3]
    std::array<int, 4> m_rotations{};


    // Размеры будущей сетки, влияют на создание таблицы m_vertices

    // Число ячеек вдоль граней (v1, v2) и (v3, v4)
    // Число ячеек вдоль граней (v1, v3) и (v2, v4)
    std::array<int, 2> m_sizes{0, 0};

    /// @brief Двумерный массив с вершинами (до вызова create_vertices
    /// может не соответствовать размерам в m_sizes).
    Array2D<BsVertex::Ptr> m_vertices{};


    // Геометрические параметры, находятся и используются при оптимизации

    /// @brief Конформный модуль блока (геометрическая характеристика)
    double m_modulus{NAN};

    /// @brief Коэффициенты для анизотропного сглаживания Лапласа
    std::array<double, 2> m_lambda{NAN, NAN};

    /// @brief Относительные размеры блока
    std::array<double, 2> m_rel_sizes{NAN, NAN};
};

} // namespace zephyr::mesh::generator
