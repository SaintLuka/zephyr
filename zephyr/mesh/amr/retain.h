/// @file Файл содержит реализацию единственной функции retain_cell, которая
/// используется для адаптации ячейки, которая не хочет адаптироваться (да, именно так).
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/faces.h>

namespace zephyr::mesh::amr {

/// @brief Функция вызывается для ячеек с флагом адаптации = 0, хотя сама
/// ячейка не разбивается и не огрубляется, у неё может измениться набор граней:
/// какие-то грани могут объединиться, а часть разбиться.
/// @param locals Хранилище локальных ячеек
/// @param ic Индекс целевой ячейки
/// @details Если какая-то грань разбивается, то adjacent у подграней
/// выставляется со старой грани. Если подграни склеиваются, то adjacent
/// у объединенной грани будет такой, какой был у старой грани.
template<int dim>
void retain_cell(AmrCells &locals, AmrCells& aliens, index_t ic) {
    // Ячейка не требует разбиения, необходимо пройти по граням,
    // возможно, необходимо разбить грань или собрать
    auto lvl_c = locals.level[ic];
    auto& adj = locals.faces.adjacent;

    for (int side = 0; side < FpC(dim); ++side) {
        auto sides = subface_sides<dim>(side);

        index_t iface = locals.face_begin[ic] + side;

        if (locals.faces.is_undefined(iface) ||
            locals.faces.is_boundary(iface)) {
            continue;
        }
#if SCRUTINY
        int rank = mpi::rank();
        if (adj.rank[iface] == rank && adj.index[iface] >= locals.size()) {
            std::cout << "AmrCell has no local neighbor through the " <<
                      side_to_string(side, dim) << " side\n";
            locals.print_info(ic);
            throw std::runtime_error("AmrCell has no local neighbor (retain_cell)");
        }
        if (adj.rank[iface] != rank &&
            (adj.alien[iface] < 0 || adj.alien[iface] >= aliens.size())) {
            std::cout << "AmrCell has no remote neighbor through the " <<
                side_to_string(side, dim) << " side\n";
            locals.print_info(ic);
            throw std::runtime_error("AmrCell has no remote neighbor (retain_cell)");
        }
#endif
        // Хранилище и индекс соседа
        auto [neibs, jc] = adj.get_neib(iface, locals, aliens);

        auto lvl_n = neibs.level[jc];

        if (lvl_c == lvl_n) {
            // Сосед имеет такой же уровень
#if SCRUTINY
            // Сложная грань вместо одинарной
            for (int i = 1; i < FpF(dim); ++i) {
                if (locals.faces.is_actual(locals.face_begin[ic] + sides[i])) {
                    std::cout << "Imbalance: lvl_c == lvl_n, complex " << side_to_string(sides[i], dim) << "\n";
                    std::cout << "MAIN CELL:\n";
                    locals.print_info(ic);
                    std::cout << "NEIB CELL:\n";
                    neibs.print_info(jc);
                    throw std::runtime_error("Imbalance: lvl_c == lvl_n, complex face");
                }
            }
#endif

            // Сосед адаптируется
            if (neibs.flag[jc] > 0) {
                split_face<dim>(ic, locals, side);
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
                if (locals.faces.is_actual(locals.face_begin[ic] + s)) {
                    ++counter;
                }
            }
            if (counter != FpF(dim)) {
                std::cout << "Current cell:\n";
                locals.print_info(ic);
                std::cout << "Neighbor:\n";
                neibs.print_info(jc);
                throw std::runtime_error("Imbalance: lvl_c < lvl_n, single face");
            }

            // Сосед не может разбиваться
            if (neibs.flag[jc] > 0) {
                std::cout << "Current cell:\n";
                locals.print_info(ic);
                std::cout << "Neighbor:\n";
                neibs.print_info(jc);
                throw std::runtime_error("Imbalance None:Split");
            }
#endif
            // Сосед огрубляется - объединяем грань
            if (neibs.flag[jc] < 0) {
                for (int s: sides) {
                    adj.rank [locals.face_begin[ic] + s] = adj.rank[iface];
                    adj.index[locals.face_begin[ic] + s] = adj.index[iface];
                    adj.alien[locals.face_begin[ic] + s] = adj.alien[iface];
                    adj.basic[locals.face_begin[ic] + s] = ic;
                }
                merge_faces<dim>(ic, locals, Side<dim>(side));
            }

            // Сосед ничего не делает - ничего не делаем
            // if (cell[adj.index].flag == 0) { }

        } else {
            // Сосед имеет уровень ниже
#if SCRUTINY
            // Сложная грань вместо одинарной
            for (int i = 1; i < FpF(dim); ++i) {

                if (locals.faces.is_actual(locals.face_begin[ic] + sides[i])) {
                    throw std::runtime_error("Imbalance: lvl_c > lvl_n, complex face");
                }
            }

            // Сосед не может огрубляться
            if (neibs.flag[jc] < 0) {
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

} // namespace zephyr::mesh::amr