// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>
#include <zephyr/mesh/amr/statistics.h>
#include <zephyr/mesh/amr/remove_undefined.h>

namespace zephyr::mesh::amr {

/// @brief Определяет положения новых ячеек (созданных при адаптации) в хранилище,
/// устанавливает параметр amrData.next.
/// @param cells Ссылка на хранилище ячеек
/// @param count Статистика адаптации
/// @details Если ячейка не изменяется, тогда next содержит индекс ячейки в
/// хранилище. Если ячейка огрубляется, тогда next содержит индекс (создаваемой)
/// родительской ячейки. Если ячейка бьется, то next содержит индекс первой
/// (создаваемой) дочерней ячейки, при этом остальные дочерние ячейки имеют
/// индексы next+1, next+2, ...
/// Алгоритм может выполняться как для всего хранилища, так и для части сетки
/// в многопроцессорном режиме. Многопоточная реализация отсутствует.
template<int dim>
void setup_positions(AmrCells &cells, const Statistics &count, const SwapLists& swap_list) {
    // TODO: Подумать над параллельной версией
    int coarse_counter = count.n_cells;
    int refine_counter = count.n_cells + count.n_parents;
    for (int ic = 0; ic < count.n_cells; ++ic) {
        if (cells.flag[ic] == 0) {
            cells.next[ic] = ic;
            continue;
        }

        if (cells.flag[ic] > 0) {
            cells.next[ic] = refine_counter;
            refine_counter += CpC(dim);
            continue;
        }

        // Главный ребенок собирает своих сиблингов
        if (cells.z_idx[ic] % CpC(dim) == 0) {
            auto sibs = get_siblings<dim>(cells, ic);
            cells.next[ic] = coarse_counter;
            for (index_t jc: sibs) {
                cells.next[jc] = coarse_counter;
            }
            ++coarse_counter;
        }
    }

    // Тождественная перестановка для новых элементов
    threads::parallel_for(
        count.n_cells, count.n_cells_large,
        [next=cells.next.data()](index_t ic) {
            next[ic] = ic;
        });

    // Наложим перестановку для новых элементов
    threads::parallel_for(
        index_t{0}, static_cast<index_t>(swap_list.actual_cells.size()),
        [next=cells.next.data(), &swap_list](index_t i) {
            next[swap_list.actual_cells[i]] = swap_list.undefined_cells[i];
        });

#if SCRUTINY
    if (coarse_counter != count.n_cells + count.n_parents || refine_counter != count.n_cells_large) {
        count.print();
        std::cout << "coarse_counter: " << coarse_counter << "\n";
        std::cout << "refine_counter: " << refine_counter << "\n";
        throw std::runtime_error("Something is wrong");
    }
#endif
}

} // namespace zephyr::mesh::amr