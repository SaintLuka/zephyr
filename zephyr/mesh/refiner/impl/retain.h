/// @file Файл содержит реализацию единственной функции retain_cell, которая
/// используется для адаптации ячейки, которая не хочет адаптироваться (да, именно так).
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/refiner/impl/common.h>
#include <zephyr/mesh/refiner/impl/faces.h>

namespace zephyr { namespace mesh { namespace impl {

using namespace ::zephyr::data;
using amrData;

/// @brief Функция вызывается для ячеек с флагом адаптации = 0, хотя сама
/// ячейка не разбивается и не огрубляется, у неё может измениться набор граней:
/// какие-то грани могут объединиться, а часть разбиться.
/// @param locals Хранилище локальных ячеек
/// @param rank Ранг текущего процесса
/// @param ic Индекс целевой ячейки
/// @details Если какая-то грань разбивается, то adjacent у подграней
/// выставляется со старой грани. Если подграни склеиваются, то adjacent
/// у объединенной грани будет такой, какой был у старой грани.
template<unsigned int dim>
void retain_cell(Storage &locals, Storage& aliens, unsigned int rank, size_t ic) {
    // Ячейка не требует разбиения, необходимо пройти по граням,
    // возможно, необходимо разбить грань или собрать
    auto cell = locals[ic];
    auto lvl_c = cell[amrData].level;

    for (int side = 0; side < FpC(dim); ++side) {
        auto sides = subface_sides<dim>(side);
        auto& face = cell[faces].list[side];
        auto adj = face.adjacent;
        if (face.is_undefined() or face.is_boundary()) {
            continue;
        }
#if SCRUTINY
        if (adj.rank == rank && adj.index >= locals.size()) {
            std::cout << "Cell has no local neighbor through the " <<
                      side_to_string(side(side % 6)) << " side\n";
            print_cell_info(cell);
                throw std::runtime_error("Cell has no local neighbor (retain_cell)");
        }
        if (adj.rank != rank && adj.ghost >= aliens.size()) {
            std::cout << "Cell has no remote neighbor through the " <<
                side_to_string(side(side % 6)) << " side\n";
            print_cell_info(cell);
            throw std::runtime_error("Cell has no remote neighbor (retain_cell)");
        }
#endif
        auto neib = adj.rank == rank ? locals[adj.index] : aliens[adj.ghost];
        auto lvl_n = neib[amrData].level;

        if (lvl_c == lvl_n) {
            // Сосед имеет такой же уровень
#if SCRUTINY
            // Сложная грань вместо одинарной
            for (int i = 1; i < FpF(dim); ++i) {
                if (cell[faces].list[sides[i]].is_actual()) {
                    throw std::runtime_error("Imbalance: lvl_c == lvl_n, complex face");
                }
            }
#endif

            // Сосед адаптируется
            if (neib[amrData].flag > 0) {
                split_face<dim>(cell, side);
            }

            // Сосед ничего не делает - ничего не делаем
            // if (neib[amrData].flag == 0) { }

            // Сосед огрубляется - ничего не делаем
            // if (neib[amrData].flag < 0) { }
        } else if (lvl_c < lvl_n) {
            // Сосед имеет уровень выше
#if SCRUTINY
            // Одинарная грань вместо сложной
            int counter = 0;
            for (int s: sides) {
                if (cell[faces].list[s].is_actual()) {
                    ++counter;
                }
            }
            if (counter != FpF(dim)) {
                throw std::runtime_error("Imbalance: lvl_c < lvl_n, single face");
            }

            // Сосед не может разбиваться
            if (neib[amrData].flag > 0) {
                std::cout << "Current cell:\n";
                print_cell_info(cell);
                std::cout << "Neighbor:\n";
                print_cell_info(neib);
                throw std::runtime_error("Imbalance None:Split");
            }
#endif
            // Сосед огрубляется - объединяем грань
            if (neib[amrData].flag < 0) {
                for (int s: sides) {
                    cell[faces].list[s].adjacent = adj;
                }
                merge_faces<dim>(cell, side(side));
            }

            // Сосед ничего не делает - ничего не делаем
            // if (cell[adj.index][amrData].flag == 0) { }

        } else {
            // Сосед имеет уровень ниже
#if SCRUTINY
            // Сложная грань вместо одинарной
            for (int i = 1; i < FpF(dim); ++i) {
                if (cell[faces].list[sides[i]].is_actual()) {
                    throw std::runtime_error("Imbalance: lvl_c > lvl_n, complex face");
                }
            }

            // Сосед не может огрубляться
            if (neib[amrData].flag < 0) {
                throw std::runtime_error("Imbalance: lvl_c > lvl_n, flag_n = coarse");
            }
#endif
            // Сосед ничего не делает - ничего не делаем
            // if (cell[adj.index][amrData].flag == 0) { }

            // Сосед адаптируется - ничего не делаем
            // if (cell[adj.index][amrData].flag == 0) { }
        }
    }
}

} // namespace impl
} // namespace mesh
} // namespace zephyr