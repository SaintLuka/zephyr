// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/siblings.h>

namespace zephyr::mesh::amr {

// Функция накладывает базовые ограничения на флаг адаптации ячейки:
//   1. Ячейка нижнего уровня не огрубляется;
//   2. Ячейка верхнего -- не разбивается;
//   3. Ячейка может иметь флаг -1 (огрубление), только если все сиблинги
//      находятся на том же процессе, имеют такой же уровень адаптации,
//      как и ячейка, а также флаг адаптации -1.
// Параметры
//   iс        -- Целевой элемент хранилища (из locals)
//   locals    -- Локальное хранилище ячеек
//   max_level -- Максимальный уровень адаптации
template<int dim>
void base_restriction(index_t ic, AmrCells &locals, int max_level) {
    scrutiny_check(ic < locals.size(), "base_restrictions: ic >= cells.size()")

    int flag = locals.flag[ic];
    // Приводим к одному из трех значений { -1, 0, 1 }
    if (flag != 0) {
        flag = flag > 0 ? 1 : -1;
    }

    int lvl = locals.level[ic];

    if (lvl + flag < 0) {
        flag = 0;
    } else {
        if (lvl + flag > max_level) {
            flag = lvl > max_level ? -1 : 0;
        }
    }

    if (flag < 0) {
        if (!can_coarse<dim>(locals, ic)) {
            flag = 0;
        }
    }

    locals.flag[ic] = flag;
}

/// @brief Выполняет функцию base_restriction для всех ячеек хранилища
/// @param max_level Максимальный уровень адаптации
template <int dim>
void base_restrictions(AmrCells &locals, int max_level) {
    threads::parallel_for(index_t{0}, index_t{locals.size()},
            base_restriction<dim>, std::ref(locals), max_level);
}

/// @brief Проверка соблюдения баланса флагов смежных ячеек и сиблингов.
/// Все сиблинги, которые хотят огрубиться, должны находиться вместе и желать
/// огрубиться тоже вместе. Уровни смежных ячеек после адаптации не должны
/// отличаться более, чем на один уровень.
/// Функция вызывается только при включенной тщательной проверке.
inline void check_flags(AmrCells& locals, AmrCells& aliens, int max_level) {
    for (index_t ic = 0; ic < locals.size(); ++ic) {
        int cell_wanted_lvl = locals.level[ic] + locals.flag[ic];

        if (cell_wanted_lvl < 0 || cell_wanted_lvl > max_level) {
            std::string message = "Wanted level (" + std::to_string(cell_wanted_lvl) + ") "
                                  + "out of range [0, " + std::to_string(max_level) + "].";
            std::cerr << message << "\n";
            throw std::runtime_error(message);
        }

        for (auto iface: locals.faces_range(ic)) {
            if (locals.faces.is_undefined(iface) ||
                locals.faces.is_boundary(iface)) {
                continue;
            }

            // Индекс соседа и хранилище соседа
            auto [neibs, jc] = locals.faces.adjacent.get_neib(iface, locals, aliens);

            int neib_wanted_lvl = neibs.level[jc] + neibs.flag[jc];
            if (std::abs(cell_wanted_lvl - neib_wanted_lvl) > 1) {
                std::string message = "Adaptation flag balance is broken.";
                std::cerr << message << "\n";
                std::cout << "cell:\n";
                locals.print_info(ic);
                std::cout << "neib:\n";
                neibs.print_info(jc);
                throw std::runtime_error(message);
            }
        }

        if (locals.flag[ic] < 0) {
            bool can = locals.dim() < 3 ?
                       can_coarse<2>(locals, ic) :
                       can_coarse<3>(locals, ic);

            if (!can) {
                std::string message = "Not all siblings want to coarse.";
                std::cerr << message << "\n";
                throw std::runtime_error(message);
            }
        }
    }
}

inline void check_flags(AmrCells& cells, int max_level) {
    AmrCells aliens;
    check_flags(cells, aliens, max_level);
}

} // namespace zephyr