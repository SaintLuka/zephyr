// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/refine.h>
#include <zephyr/mesh/amr/retain.h>
#include <zephyr/mesh/amr/coarse.h>
#include <zephyr/mesh/amr/statistics.h>

namespace zephyr::mesh::amr {

/// @brief Вызывает для ячейки соответствующий метод адаптации
/// @param locals Локальные ячейки
/// @param aliens Удаленные ячейки
/// @param rank Ранг текущего процесса
/// @param op Оператор распределения данных при огрублении и разбиении
template<int dim>
void setup_geometry_one(index_t ic, AmrCells &locals, AmrCells& aliens,
    const Distributor& op, int rank) {

    if (locals.flag[ic] == 0) {
        retain_cell<dim>(locals, aliens, ic);
        return;
    }

    if (locals.flag[ic] > 0) {
        refine_cell<dim>(locals, aliens, ic, op);
        return;
    }

    coarse_cell<dim>(locals, aliens, ic, op, rank);
}

/// @brief Осуществляет проход по ячейкам и вызывает для них
/// соответствующие методы адаптации (без MPI)
template<int dim>
void setup_geometry(AmrCells &locals, AmrCells& aliens, const Statistics &count, const Distributor& op) {
    threads::parallel_for(
            index_t{0}, index_t{count.n_cells},
            setup_geometry_one<dim>,
            std::ref(locals), std::ref(aliens), std::ref(op), 0);
}

#ifdef ZEPHYR_MPI
/// @brief Осуществляет проход по диапазону ячеек и вызывает для них
/// соответствующие методы адаптации (с MPI и без тредов)
template<int dim>
void setup_geometry(AmrCells &locals, AmrCells& aliens, const Statistics &count, const Distributor& op, int rank) {
    threads::parallel_for(
            index_t{0}, index_t{count.n_cells},
            setup_geometry_one<dim>,
            std::ref(locals), std::ref(aliens), std::ref(op), rank);
}
#endif

} // namespace zephyr::mesh::amr