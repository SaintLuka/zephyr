/// @file Файл содержит реализацию функции setup_geometry, которая создает новые ячейки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

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
void setup_geometry_one(AmrStorage::Item& cell, AmrStorage &locals, AmrStorage& aliens,
        int rank, const Distributor& op) {

    if (cell.flag == 0) {
        retain_cell<dim>(cell, locals, aliens);
        return;
    }

    if (cell.flag > 0) {
        refine_cell<dim>(cell, locals, aliens, rank, op);
        return;
    }

    coarse_cell<dim>(cell, locals, aliens, rank, op);
}

/// @brief Осуществляет проход по ячейкам и вызывает для них
/// соответствующие методы адаптации (без MPI)
template<int dim>
void setup_geometry(AmrStorage &cells, const Statistics &count, const Distributor& op) {
    static AmrStorage aliens{};

    threads::for_each<10>(
            cells.begin(), cells.begin() + count.n_cells,
            setup_geometry_one<dim>, std::ref(cells),
            std::ref(aliens), 0, std::ref(op));
}

#ifdef ZEPHYR_ENABLE_MPI
/// @brief Осуществляет проход по диапазону ячеек и вызывает для них
/// соответствующие методы адаптации (с MPI и без тредов)
template<int dim>
void setup_geometry(AmrStorage &locals, AmrStorage &aliens, int rank,
                    const Statistics &count, const Distributor& op) {
//    setup_geometry_partial<dim>(locals, aliens, rank, count, op, 0, count.n_cells);
}

#ifdef ZEPHYR_ENABLE_MULTITHREADING
/// @brief Осуществляет проход по диапазону ячеек и вызывает для них
/// соответствующие методы адаптации (с MPI и с тредами)
template<int dim>
void setup_geometry(AmrStorage &locals, AmrStorage &aliens, int rank,
        const Statistics<dim> &count, const DataDistributor& op, ThreadPool& threads) {

    auto num_tasks = threads.size();
    if (num_tasks < 2) {
        setup_geometry_partial<dim>(locals, aliens, rank, count, op, 0, count.n_cells);
        return;
    }
    std::vector<std::future<void>> results(num_tasks);

    std::int bin = count.n_cells / num_tasks + 1;
    std::int pos = 0;
    for (auto &res : results) {
        res = threads.enqueue(setup_geometry_partial<dim>,
                              std::ref(locals), std::ref(aliens), rank, std::ref(count), std::ref(op),
                              pos, std::min(pos + bin, count.n_cells)
        );
        pos += bin;
    }

    for (auto &result: results) {
        result.get();
    }
}
#endif
#endif

} // namespace zephyr::mesh::amr