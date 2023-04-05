/// @file Файл содержит реализацию функции restore_connections, которая создает связи
/// для новых ячеек, созданных после адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>
#include <zephyr/mesh/amr/statistics.h>

namespace zephyr { namespace mesh { namespace amr {

/// @brief Восстанавливает связи у части ячеек в диапазоне
/// @param locals Локальные ячейки
/// @param aliens Удаленные ячейки
/// @param count Статистика адаптации
/// @param rank Ранг процесса
/// @param from, to Диапазон ячеек (в locals)
/// @details Старые ячейки (существовшие до процедуры адаптации) в поле amrData.next
/// содержат данные о расположении ячеек, которые придут им на смену
/// (см. setup_positions). Если ячейка не меняется, тогда cells[ic].geom().next = ic.
/// Если ячейка разбивается, тогда amrData.next содержит индекс первой дочерней
/// ячейки (остальные располагаются за ней). Если ячейка огрубляется, тогда
/// amrData.next содержит индекс родительской ячейки.
/// После проведения геометрической части адаптаиции (см. setup_geometry)
/// новые ячейки уже созданы и располагаются в хранилище, но ссылаются на
/// ячейки по старым индексам. К примеру, пусть ячейка ic не меняется и имеет
/// справа соседа такого же уровня с индексом in, который хочет разбиться.
/// После проведения setup_geometry ячейка ic продолжит ссылаться на in,
/// в то время как cells[in].geom().next будет содержать индекс ребенка.
/// Таким образом, в соответствии с amrData.next ищутся актуальные ссылки на
/// реальных соседей, полученных после адаптации.
/// Грани через процессы остаются без изменений. Линковка граней между процессами
/// происходит на последнем этапе алгоритма адапатции.
template<int dim>
void restore_connections_partial(Storage &locals, Storage &aliens, int rank,
        const Statistics &count, int from, int to) {
    for (int ic = from; ic < to; ++ic) {
        auto cell = locals[ic];

        // Ячейка неактуальна, пропускаем
        if (cell.is_undefined()) continue;

        for (int i = 0; i < Faces::max_size; ++i) {
            auto &face1 = cell.geom().faces[i];
            if (face1.is_undefined() or face1.is_boundary()) continue;

            auto idx = face1.adjacent.index;
            if (face1.adjacent.rank != rank or idx > count.n_cells) {
                continue;
            }
            auto old_neib = locals[idx];

            // Актуальный сосед через грань
            if (old_neib.is_actual()) continue;

            // Сосед через грань ничего не делал - ничего не надо,
            // ссылки актуальны
            if (old_neib.flag() == 0) {
                continue;
            }

            // Сосед через грань огрубился
            if (old_neib.flag() < 1) {
                // Соседняя ячейка огрубляется
                face1.adjacent.index = old_neib.geom().next;
                continue;
            }

            // Сосед через грань адаптировался
            bool found = false;

            auto fc1 = face1.center<dim>(cell.vertices());

            for (int j = 0; j < CpC(dim); ++j) {
                auto jc = old_neib.geom().next + j;
                if (jc == ic) continue;
                if (jc > count.n_cells_large) {
                    continue;
                }

                auto neib = locals[jc];
                for (auto &face2: neib.geom().faces.list()) {
                    if (face2.is_undefined()) continue;
                    auto fc2 = face2.center<dim>(neib.vertices());

                    if ((fc1 - fc2).norm() < 1.0e-5 * cell.size()) {
                        face1.adjacent.index = jc;
                        found = true;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }
            if (!found) {
                std::cout << "Can't find neighbor through the " << side_to_string(Side(i % 6)) << " face (" << i / 6 << ")\n";
                //Storage aliens();
                //amr::print_cell_info(locals, aliens, ic);
                throw std::runtime_error("Can't find neighbor");
            }
        }
    }
}

/// @brief Устанавливает поле amrData.next в неопределенное состояние
void set_undefined_next(Storage& cells, int from, int to) {
    for (int ic = from; ic < to; ++ic) {
        cells[ic].geom().next = -1;
    }
}

/// @brief Восстанавливает связи (без MPI и без тредов)
/// @details см. restore_connections_partial
template<int dim>
void restore_connections(Storage &cells, const Statistics &count) {
    restore_connections_partial<dim>(cells, cells, 0, count, 0, count.n_cells_large);
    set_undefined_next(cells, 0, count.n_cells_large);
}

#ifdef ZEPHYR_ENABLE_MULTITHREADING
/// @brief Восстанавливает связи (без MPI с тредами)
/// @details см. restore_connections_partial
template<int dim>
void restore_connections(Storage &cells, const Statistics<dim> &count, ThreadPool& threads) {
    Storage aliens;
    auto num_tasks = threads.size();
    if (num_tasks < 2) {
        restore_connections_partial<dim>(cells, aliens, 0, count, 0, count.n_cells_large);
        set_undefined_next(cells, 0, count.n_cells_large);
        return;
    }
    std::vector<std::future<void>> results(num_tasks);

    std::int bin = count.n_cells_large / num_tasks + 1;
    std::int pos = 0;
    for (auto &res : results) {
        res = threads.enqueue(restore_connections_partial<dim>,
                              std::ref(cells), std::ref(aliens), 0, std::ref(count),
                              pos, std::min(pos + bin, count.n_cells_large)
        );
        pos += bin;
    }

    for (auto &result: results) {
        result.get();
    }

    pos = 0;
    for (auto &res : results) {
        res = threads.enqueue(set_undefined_next,
                              std::ref(cells),
                              pos, std::min(pos + bin, count.n_cells_large)
        );
        pos += bin;
    }

    for (auto &result: results) {
        result.get();
    }
}
#endif

#ifdef ZEPHYR_ENABLE_MPI
/// @brief Восстанавливает связи (с MPI и без тредов)
/// @details см. restore_connections_partial
template<int dim>
void restore_connections(Storage &locals, Storage& aliens, int rank, const Statistics<dim> &count) {
    restore_connections_partial<dim>(locals, aliens, rank, count, 0, count.n_cells_large);
    set_undefined_next(locals, 0, count.n_cells_large);
}

#ifdef ZEPHYR_ENABLE_MULTITHREADING
/// @brief Восстанавливает связи (с MPI и с тркдами)
/// @details см. restore_connections_partial
template<int dim>
void restore_connections(Storage &locals, Storage &aliens, int rank,
                         const Statistics<dim> &count, ThreadPool& threads) {
    auto num_tasks = threads.size();
    if (num_tasks < 2) {
        restore_connections_partial<dim>(locals, aliens, rank, count, 0, count.n_cells_large);
        set_undefined_next(locals, 0, count.n_cells_large);
        return;
    }
    std::vector<std::future<void>> results(num_tasks);

    std::int bin = count.n_cells_large / num_tasks + 1;
    std::int pos = 0;
    for (auto &res : results) {
        res = threads.enqueue(restore_connections_partial<dim>,
                              std::ref(locals), std::ref(aliens), rank, std::ref(count),
                              pos, std::min(pos + bin, count.n_cells_large)
        );
        pos += bin;
    }

    for (auto &result: results) {
        result.get();
    }

    pos = 0;
    for (auto &res : results) {
        res = threads.enqueue(set_undefined_next,
                              std::ref(locals),
                              pos, std::min(pos + bin, count.n_cells_large)
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