/// @file Файл содержит реализацию единственной функции retain_cell, которая
/// используется для адаптации ячейки, которая не хочет адаптироваться (да, именно так).
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr2/common.h>
#include <zephyr/mesh/amr2/faces.h>

namespace zephyr::mesh::amr2 {

/// @brief Функция вызывается для ячеек с флагом адаптации = 0, хотя сама
/// ячейка не разбивается и не огрубляется, у неё может измениться набор граней:
/// какие-то грани могут объединиться, а часть разбиться.
/// @param
/// @param cells Хранилище локальных ячеек
/// @param rank Ранг текущего процесса
/// @param ic Индекс целевой ячейки
/// @details Если какая-то грань разбивается, то adjacent у подграней
/// выставляется со старой грани. Если подграни склеиваются, то adjacent
/// у объединенной грани будет такой, какой был у старой грани.
template<int dim>
void retain_cell(index_t ic, SoaCell &cells) {
    // Ячейка не требует разбиения, необходимо пройти по граням,
    // возможно, необходимо разбить грань или собрать
    auto rank = cells.rank[ic];
    auto lvl_c = cells.level[ic];
    auto& adj = cells.faces.adjacent;

    for (int side = 0; side < FpC(dim); ++side) {
        auto sides = subface_sides<dim>(side);

        index_t iface = cells.face_begin[ic] + side;

        if (cells.faces.is_undefined(iface) or
            cells.faces.is_boundary(iface)) {
            continue;
        }
#if SCRUTINY
        if (adj.rank[iface] == rank && adj.local_index[iface] >= cells.n_locals()) {
            std::cout << "AmrCell has no local neighbor through the " <<
                      side_to_string(side % 6) << " side\n";
            cells.print_info(ic);
            throw std::runtime_error("AmrCell has no local neighbor (retain_cell)");
        }
        if (adj.rank[iface] != rank &&
            (adj.local_index[iface] < cells.n_locals() || adj.local_index[iface] >= cells.n_cells())) {
            std::cout << "AmrCell has no remote neighbor through the " <<
                side_to_string(side % 6) << " side\n";
            cells.print_info(ic);
            throw std::runtime_error("AmrCell has no remote neighbor (retain_cell)");
        }
#endif
        index_t jc = adj.local_index[iface];
        auto lvl_n = cells.level[jc];

        if (lvl_c == lvl_n) {
            // Сосед имеет такой же уровень
#if SCRUTINY
            // Сложная грань вместо одинарной
            for (int i = 1; i < FpF(dim); ++i) {
                if (cells.faces.is_actual(cells.face_begin[ic] + sides[i])) {
                    throw std::runtime_error("Imbalance: lvl_c == lvl_n, complex face");
                }
            }
#endif

            // Сосед адаптируется
            if (cells.flag[jc] > 0) {
                split_face<dim>(ic, cells, side);
            }

            // Сосед ничего не делает - ничего не делаем
            // if (neib.flag == 0) { }

            // Сосед огрубляется - ничего не делаем
            // if (neib.flag < 0) { }
        } else if (lvl_c < lvl_n) {
            // Сосед имеет уровень выше
#if SCRUTINY
            // Одинарная грань вместо сложной
            int counter = 0;
            for (int s: sides) {
                if (cells.faces.is_actual(cells.face_begin[ic] + s)) {
                    ++counter;
                }
            }
            if (counter != FpF(dim)) {
                throw std::runtime_error("Imbalance: lvl_c < lvl_n, single face");
            }

            // Сосед не может разбиваться
            if (cells.flag[jc] > 0) {
                std::cout << "Current cell:\n";
                cells.print_info(ic);
                std::cout << "Neighbor:\n";
                cells.print_info(jc);
                throw std::runtime_error("Imbalance None:Split");
            }
#endif
            // Сосед огрубляется - объединяем грань
            if (cells.flag[jc] < 0) {
                for (int s: sides) {
                    cells.faces.adjacent.rank[cells.face_begin[ic] + s] = adj.rank[iface];
                    cells.faces.adjacent.local_index[cells.face_begin[ic] + s] = adj.local_index[iface];
                    cells.faces.adjacent.owner_index[cells.face_begin[ic] + s] = adj.owner_index[iface];
                }
                merge_faces<dim>(ic, cells, Side3D(side));
            }

            // Сосед ничего не делает - ничего не делаем
            // if (cell[adj.index].flag == 0) { }

        } else {
            // Сосед имеет уровень ниже
#if SCRUTINY
            // Сложная грань вместо одинарной
            for (int i = 1; i < FpF(dim); ++i) {
                if (cells.faces.is_actual(cells.face_begin[ic] + sides[i])) {
                    throw std::runtime_error("Imbalance: lvl_c > lvl_n, complex face");
                }
            }

            // Сосед не может огрубляться
            if (cells.flag[jc] < 0) {
                throw std::runtime_error("Imbalance: lvl_c > lvl_n, flag_n = coarse");
            }
#endif
            // Сосед ничего не делает - ничего не делаем
            // if (cells.flag[jc] == 0) { }

            // Сосед адаптируется - ничего не делаем
            // if (cells.flag[jc] > 0) { }
        }
    }
}

} // namespace zephyr::mesh::amr2