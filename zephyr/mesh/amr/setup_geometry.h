/// @file Файл содержит реализацию функции setup_geometry, которая создает новые ячейки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/refine.h>
#include <zephyr/mesh/amr/retain.h>
#include <zephyr/mesh/amr/coarse.h>
#include <zephyr/mesh/amr/statistics.h>

namespace zephyr { namespace mesh { namespace amr {

/// @brief Осуществляет проход по диапазону ячеек и вызывает для них
/// соответствующие методы адаптации
/// @param locals Локальные ячейки
/// @param aliens Удаленные ячейки
/// @param rank Ранг текущего процесса
/// @param count Статистика адапатции
/// @param op Оператор распределения данных при огрублении и разбиении
/// @param from, to Диапазон ячеек
template<int dim>
void setup_geometry_partial(AmrStorage &locals, AmrStorage& aliens, int rank,
        const Statistics &count, const Distributor& op, int from, int to) {


    for (int ic = from; ic < to; ++ic) {
        AmrCell& cell = locals[ic];

        if (cell.flag == 0) {
            retain_cell<dim>(cell, locals, aliens);
            continue;
        }

        if (cell.flag > 0) {
            refine_cell<dim>(locals, aliens, rank, ic, op);
            continue;
        }

        coarse_cell<dim>(locals, aliens, rank, ic, op);
    }
}

/// @brief Осуществляет проход по диапазону ячеек и вызывает для них
/// соответствующие методы адаптации (без MPI и без тредов)
template<int dim>
void setup_geometry(AmrStorage &cells, const Statistics &count, const Distributor& op) {
    setup_geometry_partial<dim>(cells, cells, 0, count, op, 0, count.n_cells);
}

#ifdef ZEPHYR_ENABLE_MULTITHREADING
/// @brief Осуществляет проход по диапазону ячеек и вызывает для них
/// соответствующие методы адаптации (без MPI и с тредами)
template<int dim>
void setup_geometry(AmrStorage &cells, const Statistics<dim> &count,
                    const DataDistributor& op, ThreadPool& threads) {
    AmrStorage aliens;
    auto num_tasks = threads.size();
    if (num_tasks < 2) {
        setup_geometry_partial<dim>(cells, aliens, 0, count, op, 0, count.n_cells);
        return;
    }
    std::vector<std::future<void>> results(num_tasks);

    std::int bin = count.n_cells / num_tasks + 1;
    std::int pos = 0;
    for (auto &res : results) {
        res = threads.enqueue(setup_geometry_partial<dim>,
                              std::ref(cells), std::ref(aliens), 0, std::ref(count), std::ref(op),
                              pos, std::min(pos + bin, count.n_cells)
        );
        pos += bin;
    }

    for (auto &result: results) {
        result.get();
    }
}
#endif

#ifdef ZEPHYR_ENABLE_MPI
/// @brief Осуществляет проход по диапазону ячеек и вызывает для них
/// соответствующие методы адаптации (с MPI и без тредов)
template<int dim>
void setup_geometry(AmrStorage &locals, AmrStorage &aliens, int rank,
                    const Statistics<dim> &count, const DataDistributor& op) {
    setup_geometry_partial<dim>(locals, aliens, rank, count, op, 0, count.n_cells);
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

} // namespace amr
} // namespace mesh
} // namespace zephyr