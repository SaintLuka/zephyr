#pragma once

#include <memory>
#include <vector>
#include <set>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/boundary.h>

namespace zephyr::geom::generator {

class Block;
class Curve;

/// @brief Внутренний тип вершины для класса Block
struct BsVertex {
public:
    using Ptr = std::shared_ptr<BsVertex>;
    using Ref = const std::shared_ptr<BsVertex> &;

    /// @brief Конструктор
    explicit BsVertex(const Vector3d &v);

    /// @brief Конструктор
    explicit BsVertex(double x, double y);

    /// @brief Создание умного указателя
    static BsVertex::Ptr create(const Vector3d &v);

    /// @brief Создание умного указателя
    static BsVertex::Ptr create(double x, double y);

    /// @brief Число смежных вершин, больше единицы.
    /// @details Для вершины на границе области нормальным является наличие
    /// трех смежных вершин, для вершины внутри области нормальное число
    /// смежных вершин равно четырем. В этих случаях вершина считается регулярной,
    /// иначе - сингулярной. Сингулярные вершины встречаются только в углах
    /// структурированных блоков.
    int n_adjacent() const;

    /// @brief Зафиксировать вершину (удаляет смежные)
    void fix();

    /// @brief Смежные вершины
    const std::vector<BsVertex *> &adjacent_vertices() const;

    /// @brief Установить соседние вершины. Вершины массива переставляются
    /// в соответствии с правилами обхода для смежных вершин. Первая вершина
    /// остается в начале.
    void set_adjacent_vertices(const std::vector<BsVertex::Ptr> &vertices);

    /// @brief Вершина считается внутренней, если она не лежит на границе
    bool inner() const;

    /// @brief Вершина считается угловой, если она лежит на пересечении
    /// нескольких линий границ.
    bool corner() const;

    /// @brief Указатель на границу,
    /// @throw runtime_error, если более одной границы
    Curve *boundary() const;

    /// @brief Множество граничных условий, пустое множество
    /// для внутренних вершин.
    std::set<Boundary> boundaries() const;

    /// @brief Установить границу
    void add_boundary(Curve *boundary);


    int index;
    Vector3d v1;  ///< Основное положение
    Vector3d v2;  ///< Дополнительное положение

private:
    /// @brief Кривые, на которых лежит вершина.
    /// Если кривая отсутствует (nullptr), то вершина считается внутренней.
    /// Если кривая одна -- обычная вершина на границе.
    /// Если более одной -- вершина в углу.
    std::set<Curve *> m_boundaries;

    /// @brief Список смежных вершин. Используется строгий порядок заполнения.
    /// Вершины обходятся против часовой стрелки внутри области. Таким образом,
    /// для граничной вершины крайние вершины внутри списка m_adjacent также
    /// должны быть граничными.
    std::vector<BsVertex *> m_adjacent;
};

} // namespace zephyr::geom::generator
