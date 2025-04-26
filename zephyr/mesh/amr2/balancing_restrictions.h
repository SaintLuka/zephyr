/// @file Файл содержит реализацию простой, но важной функции, которая накладывает
/// базовые ограничения на флаги адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr2/common.h>
#include <zephyr/mesh/amr2/siblings.h>

namespace zephyr::mesh::amr2 {

/// @brief Функция накладывает базовые ограничения на флаг адаптации ячейки:
/// 1. Ячейка нижнего уровня не огрубляется;
/// 2. Ячейка верхнего -- не разбивается;
/// 3. Ячейка может иметь флаг -1 (огрубление), только если все сиблинги
/// находятся на том же процессе, имеют такой же уровень адаптации как и
/// ячейка, а также флаг адаптации -1.
/// @param item Целевой элемент хранилища (из locals)
/// @param cells Локальное хранилище ячеек
/// @param max_level Максимальный уровень адаптации
template<int dim>
void base_restriction(index_t ic, SoaCell &cells, int max_level) {
    scrutiny_check(ic < cells.n_locals(), "base_restrictions: ic >= cells.size()")

    int flag = cells.flag[ic];
    // Приводим к одному из трех значений { -1, 0, 1 }
    if (flag != 0) {
        flag = flag > 0 ? 1 : -1;
    }

    int lvl = cells.level[ic];

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

    cells.flag[ic] = flag;
}

/// @brief Выполняет функцию base_restriction для всех ячеек хранилища
/// @param cells Хранилище ячеек
/// @param max_level Максимальный уровень адаптации
template <int dim>
void base_restrictions(SoaCell &cells, int max_level) {
    utils::range<index_t> range(cells.n_locals());
    threads::for_each( range.begin(), range.end(),
            base_restriction<dim>, std::ref(cells), max_level);
}

/// @brief Проверка соблюдения баланса флагов смежных ячеек.
/// Уровни смежных ячеек после адаптации не должны отличаться не более,
/// чем на один уровень.
/// Функция вызывается только при включенной тщательной проверке.
void check_flags(SoaCell& cells, int max_level) {
    for (index_t ic = 0; ic < cells.n_locals(); ++ic) {
        int cell_wanted_lvl = cells.level[ic] + cells.flag[ic];

        if (cell_wanted_lvl < 0 || cell_wanted_lvl > max_level) {
            std::string message = "Wanted level (" + std::to_string(cell_wanted_lvl) + ") "
                                  + "out of range [0, " + std::to_string(max_level) + "].";
            std::cerr << message << "\n";
            throw std::runtime_error(message);
        }

        for (auto iface: cells.faces_range(ic)) {
            if (cells.faces.is_undefined(iface) ||
                cells.faces.is_boundary(iface)) {
                continue;
            }

            index_t jc = cells.faces.adjacent.local_index[iface];

            int neib_wanted_lvl = cells.level[jc] + cells.flag[jc];
            if (std::abs(cell_wanted_lvl - neib_wanted_lvl) > 1) {
                std::string message = "Adaptation flag balance is broken.";
                std::cerr << message << "\n";
                std::cout << "cell:\n";
                cells.print_info(ic);
                std::cout << "neib:\n";
                cells.print_info(jc);
                throw std::runtime_error(message);
            }
        }

        if (cells.flag[ic] < 0) {
            bool can = cells.dim < 3 ?
                       can_coarse<2>(cells, ic) :
                       can_coarse<3>(cells, ic);

            if (!can) {
                std::string message = "Not all siblings want to coarse.";
                std::cerr << message << "\n";
                throw std::runtime_error(message);
            }
        }
    }
}

} // namespace zephyr