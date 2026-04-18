#pragma once

#include <memory>
#include <vector>
#include <set>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/boundary.h>

namespace zephyr::geom::generator {

class Curve;
class BsVertex;

/// @brief Связь с соседней вершиной, если связь с внутренней вершиной,
/// тогда есть указатели на lambda1 и lambda2, если связь на границе,
/// тогда lambda2 = nullptr.
class BsEdge {
    using BsVertex_Ref = const std::shared_ptr<BsVertex>&;

    BsVertex* neib{nullptr};
    double* lambda1{nullptr};
    double* lambda2{nullptr};

public:
    /// @brief Полный конструктор
    BsEdge(BsVertex* v, double* lambda_1, double* lambda_2);

    /// @brief Создать ребро на границе
    static BsEdge Border(BsVertex_Ref v, double* lambda);

    /// @brief Создать внутреннее ребро (внутри блока)
    static BsEdge Inside(BsVertex_Ref v, double* lambda);

    /// @brief Создать ребро между парой блоков
    static BsEdge Inside(BsVertex_Ref v, double* lambda1, double* lambda2);

    /// @brief X-координата соседней вершины
    double x() const;

    /// @brief Y-координата соседней вершины
    double y() const;

    /// @brief Положение соседней вершины
    const Vector3d& pos() const;

    /// @brief Соседняя вершина на границе?
    bool boundary() const { return lambda2 == nullptr; }

    /// @brief Среднегеометрический коэффициент для сглаживания
    double lambda() const;

    /// @brief Индекс соседней вершины
    int neib_idx() const;
};


/// @brief Внутренняя вершина для класса Block
class BsVertex {
public:
    using Ptr = std::shared_ptr<BsVertex>;
    using Ref = const std::shared_ptr<BsVertex> &;

    int index;      ///< Индекс вершины
    Vector3d pos;   ///< Основное положение
    Vector3d next;  ///< Дополнительное положение

    /// @brief Конструктор
    explicit BsVertex(const Vector3d &v);

    /// @brief Конструктор
    explicit BsVertex(double x, double y);

    /// @brief Создание умного указателя
    static BsVertex::Ptr create(const Vector3d &v);

    /// @brief Создание умного указателя
    static BsVertex::Ptr create(double x, double y);

    /// @brief X-координата вершины
    double x() const { return pos.x(); }

    /// @brief Y-координата вершины
    double y() const { return pos.y(); }

    /// @brief Зафиксировать вершину (удаляет смежные)
    void fix() { m_edges.clear(); }

    /// @brief Очистить все списки (смежность и границы)
    void clear() {  m_edges.clear(); m_boundaries.clear(); }

    /// @brief Установить соседние вершины. Вершины массива переставляются
    /// в соответствии с правилами обхода для смежных вершин. Первая вершина
    /// остается в начале.
    void set_edges(const std::vector<BsEdge> &edges);

    /// @brief Число смежных вершин
    int degree() const { return static_cast<int>(m_edges.size()); }

    /// @brief Смежные вершины
    /// @details Для вершины на границе области нормальным является наличие
    /// трех смежных вершин, для вершины внутри области нормальное число
    /// смежных вершин равно четырем. В этих случаях вершина считается регулярной,
    /// иначе - сингулярной. Сингулярные вершины встречаются только в углах
    /// структурированных блоков.
    const std::vector<BsEdge> &adjacent() const { return m_edges; }

    /// @brief Вершина считается внутренней, если она не лежит на границе
    bool inner() const { return m_boundaries.empty(); }

    /// @brief Граничная вершина, если ровно одна граница
    bool is_boundary() const { return m_boundaries.size() == 1; }

    /// @brief Вершина считается угловой, если она лежит на пересечении
    /// нескольких линий границ.
    bool corner() const { return m_boundaries.size() > 1; }

    /// @brief Установить границу
    void add_boundary(Curve* boundary);

    /// @brief Указатель на границу,
    /// @throw runtime_error, если более одной границы
    Curve* boundary() const;

    /// @brief Множество граничных условий, пустое множество
    /// для внутренних вершин.
    std::set<Boundary> boundaries() const;

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
    std::vector<BsEdge> m_edges;
};

inline double BsEdge::x() const { return neib->pos.x(); }

inline double BsEdge::y() const { return neib->pos.y(); }

inline const Vector3d& BsEdge::pos() const { return neib->pos; }

inline int BsEdge::neib_idx() const { return neib->index; }

} // namespace zephyr::geom::generator
