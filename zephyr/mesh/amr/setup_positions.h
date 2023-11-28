/// @file Файл содержит реализацию функции setup_positions, которая определяет положения
/// новых ячеек (созданных при адаптации) в хранилище.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>
#include <zephyr/mesh/amr/statistics.h>

namespace zephyr { namespace mesh { namespace amr {

/// @brief Определяет положения новых ячеек (созданных при адаптации) в хранилище,
/// устанавлевает параметр amrData.next.
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
void setup_positions(AmrStorage &cells, const Statistics &count)
{
    cells.resize(count.n_cells_large);

    int coarse_counter = count.n_cells;
    int refine_counter = count.n_cells + count.n_parents;
    for (int ic = 0; ic < count.n_cells; ++ic) {
        AmrCell& cell = cells[ic];

        if (cell.flag == 0) {
            cell.next = ic;
            continue;
        }

        if (cell.flag > 0) {
            cell.next = refine_counter;
            refine_counter += CpC(dim);
            continue;
        }

        // Главный ребенок собирает своих сиблингов
        if (cell.z_idx % CpC(dim) == 0) {
            auto sibs = get_siblings<dim>(cells, ic);
            cell.next = coarse_counter;
            for (int jc: sibs) {
                cells[jc].next = coarse_counter;
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

} // namespace amr
} // namespace mesh
} // namespace zephyr