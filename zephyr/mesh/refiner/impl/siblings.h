/// @file Файл содержит набор функций для поиска сиблингов в хранилище.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <set>

#include <zephyr/mesh/refiner/impl/common.h>
#include <zephyr/mesh/refiner/impl/check.h>

namespace zephyr { namespace mesh { namespace impl {

/// @brief Сторона, по которой необходимо пройти, чтобы от одного сиблинга
/// перейти к следующему. Детали можно найти в файле _ascii.h
template <unsigned int dim>
inline std::array<side, CpC(dim)> side_to_next_sibling();

template <>
inline std::array<side, 4> side_to_next_sibling<2>() {
    return {Side::::RIGHT, Side::::TOP, Side::::BOTTOM, Side::::LEFT};
}

template <>
inline std::array<side, 8> side_to_next_sibling<3>() {
    return {
            Side::::RIGHT, Side::::TOP, Side::::FRONT, Side::::LEFT,
            Side::::BACK, Side::::LEFT, Side::::RIGHT, Side::::BOTTOM,
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
template <unsigned int dim>
bool can_coarse(Storage& cells, size_t ic) {
    std::array<side, CpC(dim)> sides = side_to_next_sibling<dim>();

    if (cells[ic][amrData].flag >= 0) {
        // Сама ячейка не хочет огрубляться
        return false;
    }

    for (unsigned int i = 0; i < CpC(dim) - 1; ++i) {
        auto cell = cells[ic];

        // локальный z-индекс
        auto z = cell[amrData].z % CpC(dim);
        auto adj = cell[faces].list[sides[z]].adjacent;

        if (adj.rank != cell[element].rank) {
            // Сосед на другом процессе
            return false;
        }

        ic = adj.index;
        auto neib = cells[ic];
        if (neib[amrData].level != cell[amrData].level) {
            // Сосед другого уровня (точно не сиблинг)
            // Может быть потомком сиблинга более высокого уровня
            return false;
        }

        if (neib[amrData].flag >= 0) {
            // Сосед не хочет огрубляться
            return false;
        }

#if SCRUTINY
        auto zc = cell[amrData].z / CpC(dim);
        auto zn = neib[amrData].z / CpC(dim);
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
template<unsigned int dim>
std::array<size_t, CpC(dim) - 1> get_siblings(Storage &cells, size_t ic) {
    std::array<size_t, CpC(dim) - 1> siblings;

    std::array<side, CpC(dim)> sides = side_to_next_sibling<dim>();

    size_t jc = ic;
    for (unsigned int i = 0; i < CpC(dim) - 1; ++i) {
        auto cell = cells[jc];

        // локальный z-индекс
        auto z = cell[amrData].z % CpC(dim);
        auto adj = cell[faces].list[sides[z]].adjacent;
        jc = adj.index;
        siblings[i] = jc;

#if SCRUTINY
        // Следующие недорозумения должны были быть устранены после выполнения
        // базовых ограничений при балансировке флагов

        // Сосед на другом процессе
        if (adj.rank != cell[element].rank) {
            throw std::runtime_error("get_siblings error: bad siblings #1");
        }

        // Сосед другого уровня
        auto neib = cells[jc];
        if (neib[amrData].level != cell[amrData].level) {
            throw std::runtime_error("get_siblings error: bad siblings #2");
        }

        auto zc = cell[amrData].z / CpC(dim);
        auto zn = neib[amrData].z / CpC(dim);
        if (zc != zn) {
            throw std::runtime_error("get_siblings error: bad siblings #3");
        }
#endif
    }

#if SCRUTINY
    /// Тестирование сиблингов

    std::set<size_t> ids;
    auto lvl = cells[ic][amrData].level;
    ids.insert(cells[ic][amrData].z % CpC(dim));

    for (size_t sib: siblings) {
        if (cells[sib][amrData].level == lvl) {
            ids.insert(cells[sib][amrData].z % CpC(dim));
        }
        else if (cells[sib][amrData].level == lvl + 1) {
            // Сиблинг через грань может уровень на единицу выше
            ids.insert((cells[sib][amrData].z / CpC(dim)) % CpC(dim));
        }
        else if (cells[sib][amrData].level == lvl + 2) {
            // Сиблинг через ребро (в 3D) или через вершину (в 2D) может иметь
            // уровень на 2 выше
            ids.insert((cells[sib][amrData].z / CpC(dim) / CpC(dim)) % CpC(dim));
        }
        else if (dim > 2 && cells[sib][amrData].level == lvl + 3) {
            // Сиблинг через вершину (в 3D) может иметь уровень на три выше
            ids.insert((cells[sib][amrData].z / CpC(dim) / CpC(dim) / CpC(dim)) % CpC(dim));
        }
        else {
            std::cout << "Current cell:\n";
            print_cell_info(cells[ic]);
            std::cout << "Sibling:\n";
            print_cell_info(cells[sib]);
            throw std::runtime_error("Different sibling levels");
        }
    }

    size_t res = 0;
    for (size_t i = 0; i < CpC(dim); ++i) {
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

} // namespace impl
} // namespace mesh
} // namespace zephyr