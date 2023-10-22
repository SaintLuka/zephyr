#pragma once

namespace zephyr::geom {

/// @brief Составной индекс смежности ячейки.
/// @details Короткое объяснение:
///               local neib           remote neib
///    rank  :    == this.rank    |    != this.rank
///    index :    < locals.size   |    < decomposition(rank).locals.size
///    ghost :    < 0             |    < aliens.size
struct Adjacent {

    /// @brief Ранг процесса, на котором находится смежная ячейка
    /// при распределенном расчете.
    int rank;

    /// @brief Индекс смежной ячейки в массиве locals (в реальном локальном
    /// хранилище или в удаленном)
    int index;

    /// @brief Индекс смежной ячейки в массиве aliens.
    int ghost;

    /// @brief Конструктор по-умолчанию.
    Adjacent() : rank(0), ghost(0), index(0) { }

    /// @brief Оператор сравнения
    bool operator!=(const Adjacent& adj) const {
        return adj.rank != rank || adj.index != index || adj.ghost != ghost;
    }
};

} // namespace zephyr::geom