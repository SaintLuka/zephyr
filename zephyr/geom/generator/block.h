#pragma once

#include <vector>
#include <memory>
#include <set>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/primitives/amr_face.h>

namespace zephyr::geom::generator {

class Curve;
class BsVertex;
class BaseVertex;

using zephyr::geom::Vector3d;
using zephyr::geom::Boundary;

using VerticesList  = std::vector<std::array<Vector3d, 4>>;
using BoundaryFlags = std::vector<std::array<Boundary, 4>>;

/// @brief Представление четырехугольного блока.
/// Основа блочно-структурированной сетки, большинство функций класса открыты
/// только для управляющего класса BlockStructured.
class Block {
private:
    using Curve_Ptr = std::shared_ptr<Curve>;
    using Curve_Ref = const std::shared_ptr<Curve> &;

    using Vertex_Ptr = std::shared_ptr<BsVertex>;
    using Vertex_Ref = const std::shared_ptr<BsVertex> &;

    using BaseVertex_Ptr = std::shared_ptr<BaseVertex>;
    using BaseVertex_Ref = const std::shared_ptr<BaseVertex> &;

public:

    /// @brief Инициализация структурированного блока.
    /// @details Вершины сортируются, чтобы получить обход против часовой
    /// стрелки, первая вершина остается на месте.
    Block &operator=(std::initializer_list<BaseVertex_Ptr> vertices);

    /// @brief Число ячеек вдоль грани
    /// @param v1, v2 Вершины грани (v1, v2)
    int size(BaseVertex_Ref v1, BaseVertex_Ref v2) const;

    /// @brief Установить число ячеек вдоль грани (v1, v2)
    void set_size(BaseVertex_Ref v1, BaseVertex_Ref v2, int N);

    /// @brief Установить границу на грань (v1, v2)
    void set_boundary(BaseVertex_Ref v1, BaseVertex_Ref v2, Curve_Ref curve);
    
    /// @brief Получить вершину
    Vertex_Ptr vertex(int i, int j) const;


private:

    /// @brief Управляющий класс
    friend class BlockStructured;

    /// @brief Пустой блок
    explicit Block(int index);

    /// @brief Индекс блока
    int index() const;

    /// @brief Разбиение вдоль граней (v1, v2) и (v3, v4)
    int size1() const;

    /// @brief Разбиение вдоль граней (v1, v4) и (v2, v3)
    int size2() const;

    /// @brief Число ячеек
    /// @param f_idx Индекс грани
    int size(int f_idx) const;

    /// @brief Смежный блок
    /// @param f_idx Индекс грани
    Block *adjacent_block(int f_idx) const;

    /// @brief Лежит ли вершина на границе
    bool is_boundary(BaseVertex_Ref v) const;

    /// @brief Граница ячейки
    /// @param f_idx Индекс грани
    Curve_Ref boundary(int f_idx) const;

    /// @brief Индекс вершины
    int vertex_index(BaseVertex_Ref v) const;

    /// @brief Базисная вершина
    /// @param v_idx Индекс вершины
    BaseVertex_Ptr base_vertex(int v_idx) const;

    /// @brief Индекс грани
    /// @param v1, v2 Вершины грани
    int face_index(BaseVertex_Ref v1, BaseVertex_Ref v2) const;

    /// @brief Индекс такой же грани у соседа
    int neib_face(int f_idx) const;

    /// @brief Связать блок с соседним
    void link(Block *block);

    /// @brief Сгенерировать узлы сетки
    void create_vertices();

    /// @brief Связать узлы сетки
    void link_vertices();

    /// @brief Сглаживание вершин в блоке
    /// @return Максимальный относительный сдивг вершин
    double smooth();

    /// @brief Обновить положение вершин
    void update();

    /// @brief Угловая вершина блока
    Vertex_Ptr &corner_vertex(int v_idx);

    /// @brief Граничная вершина блока
    Vertex_Ptr &boundary_vertex(int f_idx, int idx);

    /// @brief Вершина со второго ряда от границы
    Vertex_Ptr &preboundary_vertex(int f_idx, int idx);

    /// @brief Индекс блока
    int m_index;

    /// @brief Базовые вершины обходятся против часовой стрелки
    std::array<BaseVertex_Ptr, 4> m_base_vertices;

    /// @brief Границы области (nullptr для внутренних границ)
    std::array<Curve_Ptr, 4> m_boundaries;

    /// @brief Ссылки на соседние блоки (nullptr для границ области)
    std::array<Block *, 4> m_adjacent_blocks;

    /// @brief Поворот соседнего блока [0..3]
    std::array<int, 4> m_rotations;

    int m_size1;     ///< Число ячеек вдоль граней (v1, v2) и (v3, v4)
    int m_size2;     ///< Число ячеек вдоль вершин (v1, v4) и (v2, v3)

    /// @brief Двумерный массив с вершинами
    std::vector<std::vector<Vertex_Ptr>> m_vertices;
};

} // namespace zephyr::mesh::generator
