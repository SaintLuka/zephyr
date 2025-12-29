#pragma once

#include <vector>
#include <memory>

#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>

namespace zephyr::geom::generator {

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
    Block &operator=(std::initializer_list<BaseVertex::Ptr> vertices);

    /// @brief Число ячеек вдоль грани
    /// @param v1, v2 Вершины грани (v1, v2)
    int size(BaseVertex::Ref v1, BaseVertex::Ref v2) const;

    /// @brief Установить число ячеек вдоль грани (v1, v2)
    void set_size(BaseVertex::Ref v1, BaseVertex::Ref v2, int N);

    /// @brief Установить границу на грань (v1, v2)
    void set_boundary(BaseVertex::Ref v1, BaseVertex::Ref v2, Curve::Ref curve);
    
    /// @brief Получить вершину
    BsVertex::Ptr vertex(int i, int j) const;

    double calc_modulus() const;

    int index;

private:
    /// @brief Управляющий класс
    friend class BlockStructured;

    /// @brief Разбиение вдоль граней (v1, v2) и (v3, v4)
    int size1() const;

    /// @brief Разбиение вдоль граней (v1, v4) и (v2, v3)
    int size2() const;

    /// @brief Число ячеек
    /// @param f_idx Индекс грани
    int size(int f_idx) const;

    /// @brief Установить число ячеек вдоль грани
    void set_size_to_face(int f_idx, int N);

    /// @brief Смежный блок
    /// @param f_idx Индекс грани
    Block::Ptr adjacent_block(int f_idx) const;

    /// @brief Лежит ли вершина на границе
    bool is_boundary(BaseVertex::Ref v) const;

    /// @brief Граница ячейки
    /// @param f_idx Индекс грани
    Curve::Ref boundary(int f_idx) const;

    /// @brief Индекс вершины
    int vertex_index(BaseVertex::Ref v) const;

    /// @brief Базисная вершина
    /// @param v_idx Индекс вершины
    BaseVertex::Ptr base_vertex(int v_idx) const;

    /// @brief Индекс грани
    /// @param v1, v2 Вершины грани
    int face_index(BaseVertex::Ref v1, BaseVertex::Ref v2) const;

    /// @brief Индекс такой же грани у соседа
    int neib_face(int f_idx) const;

    /// @brief Связать блок с соседним
    void link(Block::Ref block);

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
    BsVertex::Ptr &boundary_vertex(int f_idx, int idx);

    /// @brief Вершина со второго ряда от границы
    BsVertex::Ptr &preboundary_vertex(int f_idx, int idx);


    // -------- Основные топологические параметры --------

    /// @brief Базовые вершины обходятся против часовой стрелки
    std::array<BaseVertex::Ptr, 4> m_base_vertices{}; // nullptr

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
    std::vector<std::vector<BsVertex::Ptr>> m_vertices{};
};

} // namespace zephyr::mesh::generator
