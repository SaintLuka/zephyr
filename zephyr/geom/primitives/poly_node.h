#pragma once

#include <zephyr/geom/primitives/amr_faces.h>
#include <zephyr/geom/primitives/element.h>

namespace zephyr::geom {

/// @brief Составной индекс смежной ячейки.
/// @details Короткое объяснение:
///               local neib           remote neib
///    rank  :    == this.rank    |    != this.rank
///    index :    < locals.size   |    < decomposition(rank).locals.size
///    ghost :    < 0             |    < aliens.size
struct PolyAdj {

    /// @brief Ранг процесса, на котором находится смежная ячейка
    /// при распределенном расчете.
    int rank;

    /// @brief Индекс смежной ячейки в массиве locals (в реальном локальном
    /// хранилище или в удаленном)
    int index;

    /// @brief Индекс смежной ячейки в массиве aliens.
    int ghost;

    /// @brief Индекс узла среди вершин ячейки
    short self;

    /// @brief Конструктор по-умолчанию.
    Adjacent() : rank(0), ghost(0), index(0), self_idx(0) { }

    /// @brief Оператор сравнения
    bool operator!=(const PolyAdj& adj) const {
        return adj.rank != rank || adj.index != index || adj.ghost != ghost || adj.self != self;
    }
};

/// @class Узел подвижной сетки
class PolyNode : public Element {
public:
    int max_neibs = 13; ///< Максимальное количество смежных ячеек

    Boundary boundary;  ///< Граничное условие
    Vector3d coords;    ///< Положение узла
    Vector3d shift;     ///< Смещение узла при следующем вызове move()

    /// @brief Индексы смежных ячеек
    std::array<PolyAdj, max_neibs> neibs;

    /// @brief Конструктор по умолчанию
    PolyNode() = default;

    /// @brief Переместить узел
    void move() {
        coords += shift;
        shift = {0.0, 0.0, 0.0};
    }
};

} // namespace zephyr::geom