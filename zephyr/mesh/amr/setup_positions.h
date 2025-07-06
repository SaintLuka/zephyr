/// @file Файл содержит реализацию функции setup_positions, которая определяет положения
/// новых ячеек (созданных при адаптации) в хранилище.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>
#include <zephyr/mesh/amr/statistics.h>

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
void setup_positions(AmrCells &cells, const Statistics &count)
{
    // TODO: Подумать над параллельной версией
    cells.resize_amr(count.n_cells_large);

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