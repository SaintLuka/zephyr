/// @brief Файл содержит реализацию функции setup_geometry, которая создает новые ячейки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.
#pragma once

#include <zephyr/mesh/amr2/common.h>
#include <zephyr/mesh/amr2/refine.h>
#include <zephyr/mesh/amr2/retain.h>
#include <zephyr/mesh/amr2/coarse.h>
#include <zephyr/mesh/amr2/statistics.h>

namespace zephyr::mesh::amr2 {

/// @brief Вызывает для ячейки соответствующий метод адаптации
/// @param locals Локальные ячейки
/// @param aliens Удаленные ячейки
/// @param rank Ранг текущего процесса
/// @param op Оператор распределения данных при огрублении и разбиении
template<int dim>
void setup_geometry_one(index_t ic, SoaCell &locals, SoaCell& aliens,
        int rank, const Distributor& op) {

    if (locals.flag[ic] == 0) {
        retain_cell<dim>(locals, aliens, ic);
        return;
    }

    if (locals.flag[ic] > 0) {
        refine_cell<dim>(locals, ic, rank, op);
        return;
    }

    coarse_cell<dim>(locals, aliens, ic, rank, op);
}

/// @brief Осуществляет проход по ячейкам и вызывает для них
/// соответствующие методы адаптации (без MPI)
template<int dim>
void setup_geometry(SoaCell &locals, SoaCell& aliens, const Statistics &count, const Distributor& op) {
    threads::parallel_for(index_t{0}, index_t{count.n_cells},
        setup_geometry_one<dim>,
        std::ref(locals), std::ref(aliens), 0, std::ref(op));
}

#ifdef ZEPHYR_MPI
/// @brief Осуществляет проход по диапазону ячеек и вызывает для них
/// соответствующие методы адаптации (с MPI и без тредов)
template<int dim>
void setup_geometry(AmrStorage &locals, AmrStorage &aliens, int rank,
                    const Statistics &count, const Distributor& op) {
    threads::for_each<10>(
            locals.begin(), locals.begin() + count.n_cells,
            setup_geometry_one<dim>, std::ref(locals),
            std::ref(aliens), rank, std::ref(op));
}
#endif

} // namespace zephyr::mesh::amr2