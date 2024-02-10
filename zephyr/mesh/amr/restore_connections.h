/// @file Файл содержит реализацию функции restore_connections, которая создает связи
/// для новых ячеек, созданных после адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>
#include <zephyr/mesh/amr/statistics.h>

namespace zephyr::mesh::amr {

/// @brief Восстанавливает связи у одной ячейки
/// @param cell Целевая ячейка
/// @param locals Локальные ячейки
/// @param aliens Удаленные ячейки
/// @param rank Ранг процесса
/// @param count Статистика адаптации
/// @details Старые ячейки (существовшие до процедуры адаптации) в поле amrData.next
/// содержат данные о расположении ячеек, которые придут им на смену
/// (см. setup_positions). Если ячейка не меняется, тогда cells[ic].next = ic.
/// Если ячейка разбивается, тогда amrData.next содержит индекс первой дочерней
/// ячейки (остальные располагаются за ней). Если ячейка огрубляется, тогда
/// amrData.next содержит индекс родительской ячейки.
/// После проведения геометрической части адаптаиции (см. setup_geometry)
/// новые ячейки уже созданы и располагаются в хранилище, но ссылаются на
/// ячейки по старым индексам. К примеру, пусть ячейка ic не меняется и имеет
/// справа соседа такого же уровня с индексом in, который хочет разбиться.
/// После проведения setup_geometry ячейка ic продолжит ссылаться на in,
/// в то время как cells[in].next будет содержать индекс ребенка.
/// Таким образом, в соответствии с amrData.next ищутся актуальные ссылки на
/// реальных соседей, полученных после адаптации.
/// Грани через процессы остаются без изменений. Линковка граней между процессами
/// происходит на последнем этапе алгоритма адапатции.
template<int dim>
void restore_connections_one(AmrStorage::Item& cell,
        AmrStorage &locals, AmrStorage &aliens,
        int rank, const Statistics &count) {

    // Ячейка неактуальна, пропускаем
    if (cell.is_undefined()) {
        return;
    }

    for (int i = 0; i < BFaces::max_count; ++i) {
        auto &face1 = cell.faces[i];
        if (face1.is_undefined() || face1.is_boundary()) {
            continue;
        }

        auto idx = face1.adjacent.index;
        if (face1.adjacent.rank != rank || idx > count.n_cells) {
            continue;
        }
        const auto &old_neib = locals[idx];

        // Актуальный сосед через грань
        if (old_neib.is_actual()) {
            continue;
        }

        // Сосед через грань ничего не делал - ничего не надо,
        // ссылки актуальны
        if (old_neib.flag == 0) {
            continue;
        }

        // Сосед через грань огрубился
        if (old_neib.flag < 1) {
            // Соседняя ячейка огрубляется
            face1.adjacent.index = old_neib.next;
            continue;
        }

        // Сосед через грань адаптировался
        bool found = false;

        for (int j = 0; j < CpC(dim); ++j) {
            auto jc = old_neib.next + j;

            if (jc > count.n_cells_large) {
                continue;
            }

            const auto &neib = locals[jc];
            if (&neib == &cell) {
                continue;
            }

            for (const auto &face2: neib.faces) {
                if (face2.is_undefined())
                    continue;

                if ((face1.center - face2.center).norm() < 1.0e-5 * cell.size) {
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
            std::cout << "Can't find neighbor through the " << side_to_string(Side(i % 6)) << " face (" << i / 6
                      << ")\n";
            //amr::print_cell_info(locals, aliens, ic);
            throw std::runtime_error("Can't find neighbor");
        }
    }
}

/// @brief Устанавливает поле AmrCell.next в неопределенное состояние
void set_undefined_next(AmrStorage::Item& cell) {
    cell.next = -1;
}

/// @brief Восстанавливает связи (без MPI)
/// @details см. restore_connections_partial
template<int dim>
void restore_connections(
        AmrStorage &locals, AmrStorage& aliens,
        int rank, const Statistics &count) {

    threads::for_each<20>(
            locals.begin(),
            locals.iterator(count.n_cells_large),
            restore_connections_one<dim>,
            std::ref(locals),
            std::ref(aliens),
            rank, std::ref(count)
    );

    threads::for_each(
            locals.begin(),
            locals.iterator(count.n_cells_large),
            set_undefined_next);
}

} // namespace zephyr