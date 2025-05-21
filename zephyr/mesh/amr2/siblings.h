/// @file Файл содержит набор функций для поиска сиблингов в хранилище.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <set>

#include <zephyr/mesh/amr2/common.h>

namespace zephyr::mesh::amr2 {

// Сторона, по которой необходимо пройти, чтобы от одного сиблинга перейти
// к следующему. Детали можно найти в файле _ascii.h
template <int dim>
inline constexpr std::array<Side3D, CpC(dim)> side_to_next_sibling() {
    if constexpr (dim == 2) {
        return {Side3D::RIGHT, Side3D::TOP, Side3D::BOTTOM, Side3D::LEFT};
    }
    else {
        return {
            Side3D::RIGHT, Side3D::TOP, Side3D::FRONT, Side3D::LEFT,
            Side3D::BACK, Side3D::LEFT, Side3D::RIGHT, Side3D::BOTTOM
        };
    }
}

// True если ячейка подходит для огрубления, иначе - false.
// Огрубление невозможно в следующих случаях:
//  - Хотя бы один сиблинг на другом процессе.
//  - Хотя бы один сиблинг имеет уровень выше (уже разбит).
//  - Хотя бы один сиблинг не хочет огрубляться.
// Условия, когда огрубление возможно:
//  - Все сиблинги находятся на одном процессе.
//  - Все сиблинги имеют один уровень.
//  - Все сиблинги хотят огрубиться.
template <int dim>
bool can_coarse(SoaCell& cells, int ic) {
    const auto sides = side_to_next_sibling<dim>();

    const auto& adj = cells.faces.adjacent;

    if (cells.flag[ic] >= 0) {
        // Сама ячейка не хочет огрубляться
        return false;
    }

    for (int i = 0; i < CpC(dim) - 1; ++i) {
        // локальный z-индекс
        auto z = cells.z_idx[ic] % CpC(dim);

        index_t jface = cells.face_begin[ic] + sides[z];

        if (adj.rank[jface] != cells.rank[ic]) {
            // Сосед на другом процессе
            return false;
        }

        scrutiny_check(adj.alien[jface] < 0, "can coarse, bad adjacent #1")
        scrutiny_check(adj.index[jface] >= 0, "can coarse, bad adjacent #2")
        scrutiny_check(adj.index[jface] < cells.size(), "can coarse, bad adjacent #3")

        index_t jc = adj.index[jface];
        if (cells.level[ic] != cells.level[jc]) {
            // Сосед другого уровня (точно не сиблинг)
            // Может быть потомком сиблинга более высокого уровня
            return false;
        }

        if (cells.flag[jc] >= 0) {
            // Сосед не хочет огрубляться
            return false;
        }

#if SCRUTINY
        auto zc = cells.z_idx[ic] / CpC(dim);
        auto zn = cells.z_idx[jc] / CpC(dim);
        if (zc != zn) {
            throw std::runtime_error("siblings error #1");
        }
#endif

        ic = jc;
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
std::array<int, CpC(dim) - 1> get_siblings(SoaCell &cells, index_t ic) {
    const std::array<Side3D, CpC(dim)> sides = side_to_next_sibling<dim>();

    std::array<int, CpC(dim) - 1> siblings;

    const auto& adj = cells.faces.adjacent;

    index_t jc = ic;
    for (int i = 0; i < CpC(dim) - 1; ++i) {
        // локальный z-индекс
        auto z = cells.z_idx[jc] % CpC(dim);

        index_t iface = cells.face_begin[jc] + sides[z];

#if SCRUTINY
        // Следующие недоразумения должны были быть устранены после выполнения
        // базовых ограничений при балансировке флагов

        // Сиблинг актуален и находится на другом процессе
        if (cells.is_actual(jc) && adj.rank[iface] != cells.rank[jc]) {
            throw std::runtime_error("get_siblings error: bad siblings #1");
        }

        // Сиблинг другого уровня
        if (cells.level[adj.index[iface]] != cells.level[jc]) {
            throw std::runtime_error("get_siblings error: bad siblings #2");
        }

        auto zc = cells.z_idx[jc] / CpC(dim);
        auto zn = cells.z_idx[adj.index[iface]] / CpC(dim);
        if (zc != zn) {
            throw std::runtime_error("get_siblings error: bad siblings #3");
        }
#endif

        jc = adj.index[iface];
        siblings[i] = jc;
    }

#if SCRUTINY
    /// Тестирование сиблингов

    std::set<index_t> ids;
    auto lvl = cells.level[ic];
    ids.insert(cells.z_idx[ic] % CpC(dim));

    for (index_t sib: siblings) {
        if (cells.level[sib] == lvl) {
            ids.insert(cells.z_idx[sib] % CpC(dim));
        }
        else if (cells.level[sib] == lvl + 1) {
            // Сиблинг через грань может уровень на единицу выше
            ids.insert((cells.z_idx[sib] / CpC(dim)) % CpC(dim));
        }
        else if (cells.level[sib] == lvl + 2) {
            // Сиблинг через ребро (в 3D) или через вершину (в 2D) может иметь
            // уровень на 2 выше
            ids.insert((cells.z_idx[sib] / CpC(dim) / CpC(dim)) % CpC(dim));
        }
        else if (dim > 2 && cells.level[sib] == lvl + 3) {
            // Сиблинг через вершину (в 3D) может иметь уровень на три выше
            ids.insert((cells.z_idx[sib] / CpC(dim) / CpC(dim) / CpC(dim)) % CpC(dim));
        }
        else {
            std::cout << "Current cell:\n";
            cells.print_info(ic);
            std::cout << "Sibling:\n";
            cells.print_info(sib);
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

} // namespace zephyr::mesh::amr2