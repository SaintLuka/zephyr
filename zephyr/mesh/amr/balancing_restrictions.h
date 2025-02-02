/// @file Файл содержит реализацию простой, но важной функции, которая накладывает
/// базовые ограничения на флаги адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>

namespace zephyr::mesh::amr {

/// @brief Функция накладывает базовые ограничения на флаг адаптации ячейки:
/// 1. Ячейка нижнего уровня не огрубляется;
/// 2. Ячейка верхнего -- не разбивается;
/// 3. Ячейка может иметь флаг -1 (огрубление), только если все сиблинги
/// находятся на том же процессе, имеют такой же уровень адаптации как и
/// ячейка, а также флаг адаптации -1.
/// @param item Целевой элемент хранилища (из locals)
/// @param locals Локальное хранилище ячеек
/// @param max_level Максимальный уровень адаптации
template<int dim>
void base_restriction(AmrStorage::Item& cell, AmrStorage &locals, int max_level) {
    scrutiny_check(cell.index < locals.size(), "base_restrictions: ic >= cells.size()")

    int flag = cell.flag;
    // Приводим к одному из трех значений { -1, 0, 1 }
    if (flag != 0) {
        flag = flag > 0 ? 1 : -1;
    }

    int lvl = cell.level;

    if (lvl + flag < 0) {
        flag = 0;
    } else {
        if (lvl + flag > max_level) {
            flag = lvl > max_level ? -1 : 0;
        }
    }

    if (flag < 0) {
        if (!can_coarse<dim>(cell, locals)) {
            flag = 0;
        }
    }

    cell.flag = flag;
}

/// @brief Выполняет функцию base_restriction для всех ячеек хранилища
/// @param locals Хранилище ячеек
/// @param max_level Максимальный уровень адаптации
template <int dim>
void base_restrictions(AmrStorage &locals, int max_level) {
    threads::for_each(
            locals.begin(), locals.end(),
            base_restriction<dim>, std::ref(locals), max_level);
}

/// @brief Проверка соблюдения баланса флагов смежных ячеек.
/// Уровни смежных ячеек после адаптации не должны отличаться не более,
/// чем на один уровень.
/// Функция вызывается только при включенной тщательной проверке.
void check_flags(AmrStorage& locals, AmrStorage& aliens, int max_level) {
    for(auto it = locals.begin(); it < locals.end(); ++it) {
        auto& cell = *it;

        int cell_wanted_lvl = cell.level + cell.flag;

        if (cell_wanted_lvl < 0 || cell_wanted_lvl > max_level) {
            std::string message = "Wanted level (" + std::to_string(cell_wanted_lvl) + ") "
                                  + "out of range [0, " + std::to_string(max_level) + "].";
            std::cerr << message << "\n";
            throw std::runtime_error(message);
        }

        for (auto &face: cell.faces) {
            if (face.is_undefined() || face.is_boundary()) {
                continue;
            }

            auto &adj = face.adjacent;

            auto& neib = adj.rank == cell.rank ?
                        locals[adj.index] : aliens[adj.alien];

            int neib_wanted_lvl = neib.level + neib.flag;
            if (std::abs(cell_wanted_lvl - neib_wanted_lvl) > 1) {
                std::string message = "Adaptation flag balance is broken.";
                std::cerr << message << "\n";
                throw std::runtime_error(message);
            }
        }

        if (cell.flag < 0) {
            bool can = false;
            if (cell.dim < 3) {
                can = can_coarse<2>(locals, it - locals.begin());
            }
            else{
                can = can_coarse<3>(locals, it - locals.begin());

            }
            if (!can) {
                std::string message = "Not all siblings want to coarse.";
                std::cerr << message << "\n";
                throw std::runtime_error(message);
            }
        }
    }
}

} // namespace zephyr