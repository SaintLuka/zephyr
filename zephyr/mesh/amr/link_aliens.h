/// @file Файл содержит реализацию функции link_aliens, которая использутся
/// для связывания соседей на различных процессов после выполнения функции
/// адаптации сетки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>

namespace zephyr { namespace mesh { namespace amr {

#ifdef ZEPHYR_ENABLE_MPI
/// @brief Функция вызывается перед непосредственным связыванием соседей.
/// У ячеек выставляются параметры element.rank и element.index, также
/// если ячейка имеет сосдеда на другом процессе, то выставляется
/// face.adjacent.ghost = base_id (соседа на другом процессе). После
/// пересылки соседи находятся по base_id. Функция выполняется для части
/// ячеек в хранилище.
void before_exchange_partial(
        AmrStorage& locals,
        AmrStorage& aliens,
        int rank,
        size_t from, size_t to)
{
    for (size_t ic = from; ic < to; ++ic) {
        auto cell = locals[ic];

        // Ячейка неактуальна, пропускаем
        if (cell[element].is_undefined()) continue;

        cell[element].rank = rank;
        cell[element].index = ic;

        for (auto &face1: cell.faces) {
            if (face1.is_undefined()) continue;

            if (face1.adjacent.rank == rank) {
                continue;
            }

            // Соседняя ячейка на другом процессе
            scrutiny_check(face1.adjacent.ghost < aliens.size(), "face1.adjacent.ghost > aliens.size()")

            auto neib = aliens[face1.adjacent.ghost];

            face1.adjacent.ghost = neib.b_idx;
        }
    }
}

void before_exchange(AmrStorage& locals, AmrStorage& aliens, int rank) {
    before_exchange_partial(locals, aliens, rank, 0, locals.size());
}

#ifdef ZEPHYR_ENABLE_MULTITHREADING
void before_exchange(AmrStorage& locals, AmrStorage& aliens, int rank, ThreadPool& threads) {
    auto num_tasks = threads.size();
    if (num_tasks < 2) {
        before_exchange(locals, aliens, rank);
        return;
    }
    std::vector<std::future<void>> results(num_tasks);

    std::size_t bin = locals.size() / num_tasks + 1;
    std::size_t pos = 0;
    for (auto &res : results) {
        res = threads.enqueue(before_exchange_partial,
                              std::ref(locals), std::ref(aliens), rank,
                              pos, std::min(pos + bin, locals.size())
        );
        pos += bin;
    }

    for (auto &result: results) result.get();
}
#endif

/// @brief Функция для поиска и связывания ячеек между двумя процессами.
/// Соседи ищутся по base_id, который указан в face.adjacent.ghost.
template <int dim>
void find_neighbors_partial(
        AmrStorage& locals,
        AmrStorage& aliens,
        size_t from, size_t to)
{
    for (size_t ic = from; ic < to; ++ic) {
        auto cell = locals[ic];

        // Ячейка неактуальна, пропускаем
        if (cell[element].is_undefined()) continue;

        double search_radius = 2.0 * cell[size];

        for (auto &face1: cell.faces) {
            if (face1.is_undefined() or face1.is_boundary()) continue;

            // Локальный сосед, пропускаем
            if (face1.adjacent.ghost > std::numeric_limits<int>::max()) {
                continue;
            }

            // Теперь попытаемся найти alien с таким base_id
            auto base_id = face1.adjacent.ghost;

            // Центр грани
            auto fc1 = face_center<dim>(face1, cell[vertices]);

            bool found = false;
            for (size_t in = 0; in < aliens.size(); ++in) {
                auto neib = aliens[in];
                if (neib.b_idx != base_id) continue;

                // Если neib слишком далеко, пропускаем
                if (distance(fc1, neib[coords]) > search_radius) {
                    continue;
                }

                // Обходим грани предполагаемого соседа, ищем подходящую
                for (auto &face2: neib.faces) {
                    if (face2.is_undefined()) continue;
                    auto fc2 = face_center<dim>(face2, neib[vertices]);

                    if (distance(fc1, fc2) < 1.0e-5 * cell[size]) {
                        face1.adjacent.rank = neib[element].rank;
                        face1.adjacent.index = neib[element].index;
                        face1.adjacent.ghost = in;

                        found = true;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }

            if (!found) {
                amr::print_cell_info(cell);
                std::cout << "Can't find remote neighbor\n";
                throw std::runtime_error("Can't find remote neighbor");
            }
        }
    }
}

template <int dim>
void find_neighbors(AmrStorage& locals, AmrStorage& aliens) {
    find_neighbors_partial<dim>(locals, aliens, 0, locals.size());
}

#ifdef ZEPHYR_ENABLE_MULTITHREADING
template <int dim>
void find_neighbors(AmrStorage& locals, AmrStorage& aliens, ThreadPool& threads) {
    auto num_tasks = threads.size();
    if (num_tasks < 2) {
        find_neighbors<dim>(locals, aliens);
        return;
    }
    std::vector<std::future<void>> results(num_tasks);

    std::size_t bin = locals.size() / num_tasks + 1;
    std::size_t pos = 0;
    for (auto &res : results) {
        res = threads.enqueue(find_neighbors_partial<dim>,
                              std::ref(locals), std::ref(aliens),
                              pos, std::min(pos + bin, locals.size())
        );
        pos += bin;
    }

    for (auto &result: results) result.get();
}
#endif

/// @brief Осуществляет связывание соседей на разных процессах в конце
/// выполнения функции адаптации. Подробности алгоритма можно прочитать
/// у функций before_exchange и find_neighbors
template <int dim>
void link_aliens(
        Decomposition& decomposition
        if_multithreading(, ThreadPool& threads = dummy_pool))
{
    AmrStorage& locals = decomposition.inner_elements();
    AmrStorage& aliens = decomposition.outer_elements();
    int rank = decomposition.network().rank();

    // Подготовить ячейки и face.adjacent перед обменом
    before_exchange(locals, aliens, rank if_multithreading(, threads));

    // Получили alien ячейки, они содержат данные о своем ранге и номере
    decomposition.prepare_aliens();
    decomposition.exchange();

    // Связать соседей на разных процессах
    find_neighbors<dim>(locals, aliens if_multithreading(, threads));
}

/// @brief Специализация для пустых хранилищ
template <>
void link_aliens<0>(
        Decomposition& decomposition
        if_multithreading(, ThreadPool& threads))
{
    decomposition.prepare_aliens();
    decomposition.exchange();
}
#endif

} // namespace amr
} // namespace mesh
} // namespace zephyr