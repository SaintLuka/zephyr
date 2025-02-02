/// @file Файл содержит набор функций для поиска сиблингов в хранилище.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <set>

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/check.h>

namespace zephyr { namespace mesh { namespace amr {

/// @brief Сторона, по которой необходимо пройти, чтобы от одного сиблинга
/// перейти к следующему. Детали можно найти в файле _ascii.h
template <int dim>
inline std::array<Side, CpC(dim)> side_to_next_sibling();

template <>
inline std::array<Side, 4> side_to_next_sibling<2>() {
    return {Side::RIGHT, Side::TOP, Side::BOTTOM, Side::LEFT};
}

template <>
inline std::array<Side, 8> side_to_next_sibling<3>() {
    return {
            Side::RIGHT, Side::TOP, Side::FRONT, Side::LEFT,
            Side::BACK, Side::LEFT, Side::RIGHT, Side::BOTTOM,
    };
}

/// @return True если ячейка подходит для огрубления, иначе - false.
/// @details Огрубление невозможно в следующих случаях:
///  - Хотя бы один сиблинг на другом процессе.
///  - Хотя бы один сиблинг имеет уровень выше (уже разбит).
///  - Хотя бы один сиблинг не хочет огрубляться.
/// Условия, когда огрубление возможно:
///  - Все сиблинги находятся на одном процессе.
///  - Все сиблинги имеют один уровень.
///  - Все сиблинги хотят огрубиться.
template <int dim>
bool can_coarse(AmrStorage& cells, int ic) {
    std::array<Side, CpC(dim)> sides = side_to_next_sibling<dim>();

    if (cells[ic].flag >= 0) {
        // Сама ячейка не хочет огрубляться
        return false;
    }

    for (int i = 0; i < CpC(dim) - 1; ++i) {
        auto& cell = cells[ic];

        // локальный z-индекс
        auto z = cell.z_idx % CpC(dim);
        auto adj = cell.faces[sides[z]].adjacent;

        if (adj.rank != cell.rank) {
            // Сосед на другом процессе
            return false;
        }

        ic = adj.index;
        auto& neib = cells[ic];
        if (neib.level != cell.level) {
            // Сосед другого уровня (точно не сиблинг)
            // Может быть потомком сиблинга более высокого уровня
            return false;
        }

        if (neib.flag >= 0) {
            // Сосед не хочет огрубляться
            return false;
        }

#if SCRUTINY
        auto zc = cell.z_idx / CpC(dim);
        auto zn = neib.z_idx / CpC(dim);
        if (zc != zn) {
            throw std::runtime_error("siblings error #1");
        }
#endif
    }
    return true;
}

/// @return True если ячейка подходит для огрубления, иначе - false.
/// @details Огрубление невозможно в следующих случаях:
///  - Хотя бы один сиблинг на другом процессе.
///  - Хотя бы один сиблинг имеет уровень выше (уже разбит).
///  - Хотя бы один сиблинг не хочет огрубляться.
/// Условия, когда огрубление возможно:
///  - Все сиблинги находятся на одном процессе.
///  - Все сиблинги имеют один уровень.
///  - Все сиблинги хотят огрубиться.
///
/// В фукнции используется сложная механика перехода от одного сиблинга
/// к следующему, что позволяет прервать функцию, не обходя всех сиблингов
template <int dim>
bool can_coarse(const AmrCell& main_cell, const AmrStorage& cells) {
    const std::array<Side, CpC(dim)> sides = side_to_next_sibling<dim>();

    if (main_cell.flag >= 0) {
        // Сама ячейка не хочет огрубляться
        return false;
    }

    int ic = main_cell.index;
    for (int i = 0; i < CpC(dim) - 1; ++i) {
        const auto& cell = cells[ic];

        // локальный z-индекс
        auto z = cell.z_idx % CpC(dim);
        auto adj = cell.faces[sides[z]].adjacent;

        if (adj.rank != cell.rank) {
            // Сосед на другом процессе
            return false;
        }

        ic = adj.index;
        const auto& neib = cells[ic];
        if (neib.level != cell.level) {
            // Сосед другого уровня (точно не сиблинг)
            // Может быть потомком сиблинга более высокого уровня
            return false;
        }

        if (neib.flag >= 0) {
            // Сосед не хочет огрубляться
            return false;
        }

#if SCRUTINY
        auto zc = cell.z_idx / CpC(dim);
        auto zn = neib.z_idx / CpC(dim);
        if (zc != zn) {
            throw std::runtime_error("siblings error #1");
        }
#endif
    }
    return true;
}

/// @brief Возвращает массив с индексами сиблингов, число сиблингов на единицу
/// меньше числа детей (индекс ic не добавляется в массив).
/// Предполагается, что для ячейки ic функция can_coarse возвращает true,
/// в обратном случае поведение функции неоопределено.
/// @param cells Хранилище ячеек
/// @param ic Целевая ячейка (от которой запрос)
template<int dim>
std::array<int, CpC(dim) - 1> get_siblings(AmrStorage &cells, int ic) {
    const std::array<Side, CpC(dim)> sides = side_to_next_sibling<dim>();

    std::array<int, CpC(dim) - 1> siblings;

    int jc = ic;
    for (int i = 0; i < CpC(dim) - 1; ++i) {
        auto& cell = cells[jc];

        // локальный z-индекс
        auto z = cell.z_idx % CpC(dim);
        auto adj = cell.faces[sides[z]].adjacent;
        jc = adj.index;
        siblings[i] = jc;

#if SCRUTINY
        // Следующие недорозумения должны были быть устранены после выполнения
        // базовых ограничений при балансировке флагов

        // Сиблинг актуален и находится на другом процессе
        if (cell.is_actual() && adj.rank != cell.rank) {
            throw std::runtime_error("get_siblings error: bad siblings #1");
        }

        // Сиблинг другого уровня
        auto& neib = cells[jc];
        if (neib.level != cell.level) {
            throw std::runtime_error("get_siblings error: bad siblings #2");
        }

        auto zc = cell.z_idx / CpC(dim);
        auto zn = neib.z_idx / CpC(dim);
        if (zc != zn) {
            throw std::runtime_error("get_siblings error: bad siblings #3");
        }
#endif
    }

#if SCRUTINY
    /// Тестирование сиблингов

    std::set<int> ids;
    auto lvl = cells[ic].level;
    ids.insert(cells[ic].z_idx % CpC(dim));

    for (int sib: siblings) {
        if (cells[sib].level == lvl) {
            ids.insert(cells[sib].z_idx % CpC(dim));
        }
        else if (cells[sib].level == lvl + 1) {
            // Сиблинг через грань может уровень на единицу выше
            ids.insert((cells[sib].z_idx / CpC(dim)) % CpC(dim));
        }
        else if (cells[sib].level == lvl + 2) {
            // Сиблинг через ребро (в 3D) или через вершину (в 2D) может иметь
            // уровень на 2 выше
            ids.insert((cells[sib].z_idx / CpC(dim) / CpC(dim)) % CpC(dim));
        }
        else if (dim > 2 && cells[sib].level == lvl + 3) {
            // Сиблинг через вершину (в 3D) может иметь уровень на три выше
            ids.insert((cells[sib].z_idx / CpC(dim) / CpC(dim) / CpC(dim)) % CpC(dim));
        }
        else {
            std::cout << "Current cell:\n";
            cells[ic].print_info();
            std::cout << "Sibling:\n";
            cells[sib].print_info();
            throw std::runtime_error("Different sibling levels");
        }
    }

    int res = 0;
    for (int i = 0; i < CpC(dim); ++i) {
        res += ids.count(i);
    }
    if (res != CpC(dim)) {
        std::cout << "z_loc: ";
        for (auto kek: ids) {
            std::cout << kek << " ";
        }
        std::cout << "\n";
        throw std::runtime_error("Strange siblings set");
    }
#endif

    return siblings;
}

/// @brief Возвращает массив с индексами сиблингов, число сиблингов на единицу
/// меньше числа детей (индекс самой ячейки не добавляется).
/// Предполагается, что для главной ячейки функция can_coarse возвращает true,
/// в обратном случае поведение функции неоопределено.
/// @param main_cell Целевая ячейка (от которой запрос)
/// @param locals Хранилище ячеек
template<int dim>
std::array<int, CpC(dim) - 1> get_siblings(const AmrCell& main_cell, AmrStorage &locals) {
    const std::array<Side, CpC(dim)> sides = side_to_next_sibling<dim>();

    std::array<int, CpC(dim) - 1> siblings;

    int jc = main_cell.index;
    for (int i = 0; i < CpC(dim) - 1; ++i) {
        const auto& cell = locals[jc];

        // локальный z-индекс
        auto z = cell.z_idx % CpC(dim);
        jc = cell.faces[sides[z]].adjacent.index;
        siblings[i] = jc;

#if SCRUTINY
        // Следующие недорозумения должны были быть устранены после выполнения
        // базовых ограничений при балансировке флагов
        auto adj = cell.faces[sides[z]].adjacent;

        // Сиблинг актуален и находится на другом процессе
        if (cell.is_actual() && adj.rank != cell.rank) {
            throw std::runtime_error("get_siblings error: bad siblings #1");
        }

        // Сиблинг другого уровня
        auto& neib = locals[jc];
        if (neib.level != cell.level) {
            throw std::runtime_error("get_siblings error: bad siblings #2");
        }

        auto zc = cell.z_idx / CpC(dim);
        auto zn = neib.z_idx / CpC(dim);
        if (zc != zn) {
            throw std::runtime_error("get_siblings error: bad siblings #3");
        }
#endif
    }

#if SCRUTINY
    /// Тестирование сиблингов

    std::set<int> ids;
    auto lvl = main_cell.level;
    ids.insert(main_cell.z_idx % CpC(dim));

    for (int sib: siblings) {
        if (locals[sib].level == lvl) {
            ids.insert(locals[sib].z_idx % CpC(dim));
        }
        else if (locals[sib].level == lvl + 1) {
            // Сиблинг через грань может уровень на единицу выше
            ids.insert((locals[sib].z_idx / CpC(dim)) % CpC(dim));
        }
        else if (locals[sib].level == lvl + 2) {
            // Сиблинг через ребро (в 3D) или через вершину (в 2D) может иметь
            // уровень на 2 выше
            ids.insert((locals[sib].z_idx / CpC(dim) / CpC(dim)) % CpC(dim));
        }
        else if (dim > 2 && locals[sib].level == lvl + 3) {
            // Сиблинг через вершину (в 3D) может иметь уровень на три выше
            ids.insert((locals[sib].z_idx / CpC(dim) / CpC(dim) / CpC(dim)) % CpC(dim));
        }
        else {
            std::cout << "Current cell:\n";
            main_cell.print_info();
            std::cout << "Sibling:\n";
            locals[sib].print_info();
            throw std::runtime_error("Different sibling levels");
        }
    }

    int res = 0;
    for (int i = 0; i < CpC(dim); ++i) {
        res += ids.count(i);
    }
    if (res != CpC(dim)) {
        std::cout << "z_loc: ";
        for (auto kek: ids) {
            std::cout << kek << " ";
        }
        std::cout << "\n";
        throw std::runtime_error("Strange siblings set");
    }
#endif

    return siblings;
}

} // namespace amr
} // namespace mesh
} // namespace zephyr