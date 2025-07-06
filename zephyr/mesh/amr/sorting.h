/// @file Реализация сортировки ячеек в хранилище
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>

namespace zephyr::mesh::amr {

/// @brief Функция сортировки ячеек с сохранением связности
/// @details Не самый быстрый алгоритм, нужен для проверки гипотез
/// по поводу производительности других алгоритмов
void sorting(AmrCells& cells) {
    throw std::runtime_error("NOt implemented sorting");

    /*
    size_t n_cells = cells.size();

    std::vector<size_t> map1(n_cells);
    for (size_t ic = 0; ic < n_cells; ++ic) {
        map1[ic] = ic;
    }
    std::sort(map1.begin(), map1.end(), [&cells](size_t i, size_t j) -> bool {
        auto& cell_i = cells[i];
        auto& cell_j = cells[j];

        if (cell_i.level != cell_j.level) {
            return cell_i.level < cell_j.level;
        }

        if (cell_i.b_idx != cell_j.b_idx) {
            return cell_i.b_idx < cell_j.b_idx;
        }

        return cell_i.z_idx < cell_j.z_idx;
    });

    std::vector<size_t> map(n_cells);
    for (size_t ic = 0; ic < n_cells; ++ic) {
        map[map1[ic]] = ic;
    }

    for (size_t ic = 0; ic < n_cells; ++ic) {
        for (auto &face: cells[ic].faces) {
            if (face.is_undefined()) continue;
            if (face.adjacent.index > n_cells) continue;

            face.adjacent.index = map[face.adjacent.index];
        }
    }

    size_t ic = 0;
    while(ic < n_cells) {
        size_t jc = map[ic];
        if (ic != jc) {
            throw std::runtime_error("ZAZAZAZA");

            _item_ item_i = cells[ic][item];
            _item_ item_j = cells[jc][item];

            std::vector<char*> temp(item_j.size());
            memcpy(temp.data(), item_j.pos(), item_i.size());
            memcpy(item_j.pos(), item_i.pos(), item_i.size());
            memcpy(item_i.pos(), temp.data(), item_i.size());


            //std::swap(map[ic], map[jc]);
            map[ic] = map[jc];
            map[jc] = jc;
        }
        else {
            ++ic;
        }
    }
    */
}

} // namespace zephyr