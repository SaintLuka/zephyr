/// @file Файл содержит реализацию простой, но важной функции, которая накладывает
/// базовые ограничения на флаги адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>

#ifdef ZEPHYR_ENABLE_MULTITHREADING
#include <zephyr/multithreading/thread-pool.h>
#endif

namespace zephyr { namespace mesh { namespace amr {

/// @brief Функция осуществляет обход по части ячеек хранилища и накладывает
/// базовые ограничения на флаги адаптации:
/// 1. Ячейка нижнего уровня не огрубляется;
/// 2. Ячейка верхнего -- не разбивается;
/// 3. Ячейка может иметь флаг -1 (огрубление), только если все сиблинги
/// находятся на том же процессе, имеют такой же уровень адаптации как и
/// ячейка, а также флаг адаптации -1.
/// @param cells Хранилище ячеек
/// @param max_level Максимальный уровень адаптации
/// @param from, to Диапазон ячеек для обхода
template<int dim>
void base_restrictions_partial(Storage &cells, int max_level, int from, int to) {
    for (int ic = from; ic < to; ++ic) {
        scrutiny_check(ic < cells.size(), "base_restrictions: ic >= cells.size()")
        auto cell = cells[ic];

        short flag = cell.flag();
        // Приводим к одному из трех значений { -1, 0, 1 }
        if (flag != 0) {
            flag = flag > 0 ? 1 : -1;
        }

        int lvl = cell.level();

        if (lvl + flag < 0) {
            flag = 0;
        } else {
            if (lvl + flag > max_level) {
                flag = lvl > max_level ? -1 : 0;
            }
        }

        if (flag < 0) {
            if (!can_coarse<dim>(cells, ic)) {
                flag = 0;
            }
        }

        cell.geom().flag = flag;
    }
}

/// @brief Выполняет функцию base_restriction_partial для всех ячеек хранилища
/// в однопоточном режиме
/// @param cells Хранилище ячеек
/// @param max_level Максимальный уровень адаптации
template <int dim>
void base_restrictions(Storage &cells, int max_level) {
    base_restrictions_partial<dim>(cells, max_level, 0, cells.size());
}

#ifdef ZEPHYR_ENABLE_MULTITHREADING
/// @brief Выполняет функцию base_restriction_partial для всех ячеек хранилища
/// в многопоточном режиме
/// @param cells Хранилище ячеек
/// @param max_level Максимальный уровень адаптации
/// @param threads Ссылка ну пул тредов
template <int dim>
void base_restrictions(Storage &cells, int max_level, ThreadPool& threads) {
    auto num_tasks = threads.size();
    if (num_tasks < 2) {
        // Вызов однопоточной версии
        base_restrictions_partial<dim>(cells, max_level, 0, cells.size());
        return;
    }
    std::vector<std::future<void>> results(num_tasks);

    std::int bin = cells.size() / num_tasks + 1;
    std::int pos = 0;
    for (auto &res : results) {
        res = threads.enqueue(base_restrictions_partial<dim>,
                              std::ref(cells), max_level,
                              pos, std::min(pos + bin, cells.size())
        );
        pos += bin;
    }
    for (auto &result: results)
        result.get();
}
#endif

void check_flags(Storage& locals, int max_level, Storage& aliens) {
    for(auto cell: locals) {
        int cell_wanted_lvl = cell.level() + cell.flag();

        if (cell_wanted_lvl < 0 || cell_wanted_lvl > max_level) {
            std::string message = "Wanted level out of range [0, " + std::to_string(max_level) + "].";
            std::cerr << message << "\n";
            throw std::runtime_error(message);
        }

        for (auto &face: cell.geom().faces.list()) {
            if (face.is_undefined() or face.is_boundary()) {
                continue;
            }

            auto &adj = face.adjacent;

            auto neib = adj.rank == cell.rank() ?
                        locals[adj.index] : aliens[adj.ghost];

            int neib_wanted_lvl = neib.level() + neib.flag();
            if (std::abs(cell_wanted_lvl - neib_wanted_lvl) > 1) {
                std::string message = "Adaptation flag balance is broken.";
                std::cerr << message << "\n";
                throw std::runtime_error(message);
            }
        }

        if (cell.flag() < 0) {
            bool can = false;
            if (cell.dim() < 3) {
                can = can_coarse<2>(locals, cell - locals.begin());
            }
            else{
                can = can_coarse<3>(locals, cell - locals.begin());

            }
            if (!can) {
                std::string message = "Not all siblings want to coarse.";
                std::cerr << message << "\n";
                throw std::runtime_error(message);
            }
        }
    }
}

} // namespace amr
} // namespace mesh
} // namespace zephyr