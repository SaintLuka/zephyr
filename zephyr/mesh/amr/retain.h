// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/faces.h>

namespace zephyr::mesh::amr {

/// @brief Функция вызывается для ячеек с флагом адаптации = 0, хотя сама ячейка
/// не разбивается и не огрубляется, у неё может измениться набор граней:
/// какие-то грани могут объединиться, а какие-то разбиться.
/// @param locals Хранилище локальных ячеек
/// @param aliens Хранилище ячеек с других процессов
/// @param ic Индекс целевой ячейки в locals
template<int dim>
void retain_cell(AmrCells &locals, AmrCells& aliens, index_t ic) {
    // Ячейка не требует разбиения, необходимо пройти по граням,
    // возможно, необходимо разбить грань или собрать
    auto lvl_c = locals.level[ic];
    auto& adj = locals.faces.adjacent;

    for (auto side: Side<dim>::items()) {
        index_t face_beg = locals.face_begin[ic];
        index_t iface = face_beg + side;

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
            auto subfaces = side.subfaces();
            // Сложная грань вместо одинарной
            for (int i = 1; i < FpF(dim); ++i) {
                if (locals.faces.is_actual(face_beg + subfaces[i])) {
                    std::cout << "Imbalance: lvl_c == lvl_n, complex " << side_to_string(subfaces[i], dim) << "\n";
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
                split_face<dim>(locals, ic, side);

                index_t neib_next = neibs.next[jc];
                for (auto subface: side.subfaces()) {
                    // TODO: MPI VERSION
                    adj.index[face_beg + subface] = neib_next + subface.neib_child();
                }
            }

            // Сосед ничего не делает - ничего не делаем
            // if (neib.flag == 0) { }

            // Сосед огрубляется
            if (neibs.flag[jc] < 0) {
                // TODO: MPI VERSION
                adj.index[iface] = neibs.next[jc];
            }
        }
        else if (lvl_c < lvl_n) {
            // Сосед имеет уровень выше
#if SCRUTINY
            // Одинарная грань вместо сложной
            int counter = 0;
            for (int subface: side.subfaces()) {
                if (locals.faces.is_actual(face_beg + subface)) {
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
                merge_faces<dim>(locals, ic, side);
                // TODO: MPI VERSION
                adj.index[iface] = neibs.next[jc];
            }

            // Сосед ничего не делает - ничего не делаем
            // if (cell[adj.index].flag == 0) { }

        } else {
            // Сосед имеет уровень ниже
#if SCRUTINY
            auto sides = side.subfaces();
            // Сложная грань вместо одинарной
            for (int i = 1; i < FpF(dim); ++i) {
                if (locals.faces.is_actual(face_beg + sides[i])) {
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

            // Сосед адаптируется
            if (neibs.flag[jc] > 0) {
                // TODO: MPI VERSION
                adj.index[iface] = neibs.next[jc] + side.adjacent_child(locals.z_idx[ic] % CpC(dim));
            }
        }
    }
}

} // namespace zephyr::mesh::amr