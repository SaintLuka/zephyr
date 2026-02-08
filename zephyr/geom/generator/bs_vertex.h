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
class BsVertex {
public:
    using Ptr = std::shared_ptr<BsVertex>;
    using Ref = const std::shared_ptr<BsVertex> &;

    /// @brief Связь с соседней вершиной, если связь с внутренней вершиной,
    /// тогда есть указатели на ratio1 и ratio2, если связь на границе,
    /// тогда ratio2 = nullptr.
    struct Edge {
        BsVertex* neib{nullptr};
        double* ratio1{nullptr};
        double* ratio2{nullptr};

        bool boundary() const {
            return ratio2 == nullptr;
        }

        double ratio() const {
            z_assert(ratio1 != nullptr, "ratio1 is null");
            if (boundary()) { return (*ratio1); }
            return std::sqrt((*ratio1) * (*ratio2));
        }
    };

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
    int degree() const;

    /// @brief Зафиксировать вершину (удаляет смежные)
    void fix();

    /// @brief Очистить все списки
    void clear();

    /// @brief Смежные вершины
    const std::vector<Edge> &adjacent() const;

    /// @brief Установить соседние вершины. Вершины массива переставляются
    /// в соответствии с правилами обхода для смежных вершин. Первая вершина
    /// остается в начале.
    void set_edges(const std::vector<Edge> &edges);

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

    double x() const { return v1.x(); }
    double y() const { return v1.y(); }


    int index;
    Vector3d v1;  ///< Основное положение
    Vector3d v2;  ///< Дополнительное положение

private:
    /// @brief Кривые, на которых лежит вершина.
    /// Если кривая отсутствует (nullptr), то вершина считается внутренней.
    /// Если кривая одна -- обычная вершина на границе.
    /// Если более одной -- вершина в углу, должна быть неподвижна.
    std::set<Curve*> m_boundaries;

    /// @brief Список смежных вершин. Используется строгий порядок заполнения.
    /// Вершины обходятся против часовой стрелки внутри области. Таким образом,
    /// для граничной вершины крайние вершины внутри списка m_adjacent также
    /// должны быть граничными.
    std::vector<Edge> m_edges;
};

} // namespace zephyr::geom::generator
