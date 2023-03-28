/// @file Файл содержит реализацию функции setup_positions, которая определяет положения
/// новых ячеек (созданных при адаптации) в хранилище.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/refiner/impl/common.h>
#include <zephyr/mesh/refiner/impl/siblings.h>
#include <zephyr/mesh/refiner/impl/statistics.h>

namespace zephyr { namespace mesh { namespace impl {

using namespace ::zephyr::data;
using amrData;

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
template<unsigned int dim>
void setup_positions(Storage &cells, const Statistics<dim> &count
    if_multithreading(, ThreadPool& threads = dummy_pool))
{
    cells.resize(count.n_cells_large);

    size_t coarse_counter = count.n_cells;
    size_t refine_counter = count.n_cells + count.n_parents;
    for (size_t ic = 0; ic < count.n_cells; ++ic) {
        auto cell = cells[ic];

        if (cell[amrData].flag == 0) {
            cell[amrData].next = ic;
            continue;
        }

        if (cell[amrData].flag > 0) {
            cell[amrData].next = refine_counter;
            refine_counter += CpC(dim);
            continue;
        }

        // Главный ребенок собирает своих сиблингов
        if (cells[ic][amrData].z % CpC(dim) == 0) {
            auto sibs = get_siblings<dim>(cells, ic);
            cells[ic][amrData].next = coarse_counter;
            for (size_t jc: sibs) {
                cells[jc][amrData].next = coarse_counter;
            }
            ++coarse_counter;
        }
    }

    if (coarse_counter != count.n_cells + count.n_parents || refine_counter != count.n_cells_large) {
        count.print();
        std::cout << "coarse_counter: " << coarse_counter << "\n";
        std::cout << "refine_counter: " << refine_counter << "\n";
        throw std::runtime_error("Something is wrong");
    }
}

} // namespace impl
} // namespace mesh
} // namespace zephyr