#pragma once

namespace zephyr::geom {

/// @brief Составной индекс смежности.
/// @details Короткое объяснение:
///               local neib           remote neib
///    rank  :    == this.rank    |    != this.rank
///    index :    < locals.size   |    < decomposition(rank).locals.size
///    alien :    < 0             |    < aliens.size
struct Adjacent {

    /// @brief Ранг процесса, на котором находится смежная ячейка
    /// при распределенном расчете.
    int rank = -1;

    /// @brief Индекс смежной ячейки в массиве locals (в реальном локальном
    /// хранилище или в удаленном)
    int index = -1;

    /// @brief Индекс смежной ячейки в массиве aliens.
    int alien = -1;

    /// @brief Конструктор по-умолчанию.
    Adjacent() = default;

    /// @brief Оператор сравнения
    bool operator!=(const Adjacent& adj) const {
        return adj.rank != rank || adj.index != index || adj.alien != alien;
    }
};

} // namespace zephyr::geom