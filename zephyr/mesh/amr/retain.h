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

        if (locals.faces.is_undefined(iface)) {
            continue;
        }
        if (locals.faces.is_boundary(iface)) {
            // Сама ячейка может переезжать
            scrutiny_check(locals.simple_face(ic, side), "Complex boundary face");
            adj.index[iface] = locals.next[ic];
            adj.basic[iface] = locals.next[ic];
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
        // Хранилище и индекс соседа, если соседи более высокого уровня,
        // то это какой-то из соседей
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
            index_t neib_next = neibs.next[jc];

            if (neibs.flag[jc] > 0) {
                // Сосед адаптируется
                split_face<dim>(locals, ic, side);

                for (auto subface: side.subfaces()) {
                    // TODO: MPI VERSION
                    adj.index[face_beg + subface] = locals.next[neib_next + subface.neib_child()];
                    adj.basic[face_beg + subface] = locals.next[ic];
                }
            }
            else {
                adj.basic[iface] = locals.next[ic];

                if (neibs.flag[jc] == 0) {
                    // Сосед ничего не делает, но может переехать
                    // TODO: MPI VERSION
                    adj.index[iface] = neib_next;
                }
                else {
                    // Сосед огрубляется (neib_flag < 0)
                    // TODO: MPI VERSION
                    adj.index[iface] = locals.next[neib_next];
                }
            }
        }
        else if (lvl_c < lvl_n) {
            // Соседи (их два/четыре) имеют уровень выше
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
            if (neibs.flag[jc] < 0) {
                // Соседи огрубляются - объединяем грань
                merge_faces<dim>(locals, ic, side);

                // TODO: MPI VERSION
                adj.index[iface] = locals.next[neibs.next[jc]];
                adj.basic[iface] = locals.next[ic];
            }
            else {
                // Сосед ничего не делает, но может переехать
                for (auto subface: side.subfaces()) {
                    index_t jface = face_beg + subface;
                    // Хранилище и индекс соседа
                    auto [neibs2, kc] = adj.get_neib(jface, locals, aliens);

                    // TODO: MPI VERSION
                    adj.index[jface] = neibs2.next[kc];
                    adj.basic[jface] = locals.next[ic];
                }
            }
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
            if (neibs.flag[jc] == 0) {
                // Сосед ничего не делает, но может переехать
                // TODO: MPI VERSION
                adj.index[iface] = neibs.next[jc];
                adj.basic[iface] = locals.next[ic];
            }
            else {
                // Сосед адаптируется (neib_flag > 0)
                // TODO: MPI VERSION
                index_t neib_next = neibs.next[jc] + side.adjacent_child(locals.z_idx[ic] % CpC(dim));
                adj.index[iface] = locals.next[neib_next];
                adj.basic[iface] = locals.next[ic];
            }
        }
    }
}

} // namespace zephyr::mesh::amr