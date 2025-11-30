// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/faces.h>

namespace zephyr::mesh::amr {

/// @brief Добавить связь, (i, j) - координаты дочерней ячейки, (x, y) - вектор направления
template<int i, int j, int x, int y>
void make_link_2D(AmrCells& cells, index_t main_child, int rank) {
    index_t iface = cells.face_begin[main_child + Quad::iss<i, j>()] + Side2D::by_dir<x, y>();
    cells.faces.boundary[iface] = Boundary::ORDINARY;
    cells.faces.adjacent.rank[iface] = rank;
    cells.faces.adjacent.index[iface] = cells.next[main_child + Quad::iss<i + 2 * x, j + 2 * y>()];
    cells.faces.adjacent.alien[iface] = -1;
}

/// @brief Добавить связь, (i, j, k) - координаты дочерней ячейки, (x, y, z) - вектор направления
template<int i, int j, int k, int x, int y, int z>
void make_link_3D(AmrCells& cells, index_t main_child, int rank) {
    index_t iface = cells.face_begin[main_child + Cube::iss<i, j, k>()] +  + Side3D::by_dir<x, y, z>();
    cells.faces.boundary[iface] = Boundary::ORDINARY;
    cells.faces.adjacent.rank[iface] = rank;
    cells.faces.adjacent.index[iface] = cells.next[main_child + Cube::iss<i + 2 * x, j + 2 * y, k + 2 * z>()];
    cells.faces.adjacent.alien[iface] = -1;
}

/// @brief Связать грани соседних дочерних ячеек (внутри родительской). Помним,
/// что дочерние ячейки располагаются последовательно, нумерация начнется с ic.
/// @param cells Локальное хранилище ячеек
/// @param ic Индекс первой (главной) ячейки
template <int dim>
void link_siblings(AmrCells& cells, index_t ic, int rank) {
    if constexpr (dim == 2) {
        make_link_2D<-1, -1, +1, 0>(cells, ic, rank);
        make_link_2D<-1, -1, 0, +1>(cells, ic, rank);

        make_link_2D<+1, -1, -1, 0>(cells, ic, rank);
        make_link_2D<+1, -1, 0, +1>(cells, ic, rank);

        make_link_2D<-1, +1, +1, 0>(cells, ic, rank);
        make_link_2D<-1, +1, 0, -1>(cells, ic, rank);

        make_link_2D<+1, +1, -1, 0>(cells, ic, rank);
        make_link_2D<+1, +1, 0, -1>(cells, ic, rank);
    }
    else {
        make_link_3D<-1, -1, -1, +1, 0, 0>(cells, ic, rank);
        make_link_3D<-1, -1, -1, 0, +1, 0>(cells, ic, rank);
        make_link_3D<-1, -1, -1, 0, 0, +1>(cells, ic, rank);

        make_link_3D<+1, -1, -1, -1, 0, 0>(cells, ic, rank);
        make_link_3D<+1, -1, -1, 0, +1, 0>(cells, ic, rank);
        make_link_3D<+1, -1, -1, 0, 0, +1>(cells, ic, rank);

        make_link_3D<-1, +1, -1, +1, 0, 0>(cells, ic, rank);
        make_link_3D<-1, +1, -1, 0, -1, 0>(cells, ic, rank);
        make_link_3D<-1, +1, -1, 0, 0, +1>(cells, ic, rank);

        make_link_3D<+1, +1, -1, -1, 0, 0>(cells, ic, rank);
        make_link_3D<+1, +1, -1, 0, -1, 0>(cells, ic, rank);
        make_link_3D<+1, +1, -1, 0, 0, +1>(cells, ic, rank);

        make_link_3D<-1, -1, +1, +1, 0, 0>(cells, ic, rank);
        make_link_3D<-1, -1, +1, 0, +1, 0>(cells, ic, rank);
        make_link_3D<-1, -1, +1, 0, 0, -1>(cells, ic, rank);

        make_link_3D<+1, -1, +1, -1, 0, 0>(cells, ic, rank);
        make_link_3D<+1, -1, +1, 0, +1, 0>(cells, ic, rank);
        make_link_3D<+1, -1, +1, 0, 0, -1>(cells, ic, rank);

        make_link_3D<-1, +1, +1, +1, 0, 0>(cells, ic, rank);
        make_link_3D<-1, +1, +1, 0, -1, 0>(cells, ic, rank);
        make_link_3D<-1, +1, +1, 0, 0, -1>(cells, ic, rank);

        make_link_3D<+1, +1, +1, -1, 0, 0>(cells, ic, rank);
        make_link_3D<+1, +1, +1, 0, -1, 0>(cells, ic, rank);
        make_link_3D<+1, +1, +1, 0, 0, -1>(cells, ic, rank);
    }
}

/// @brief Проверить связи между дочерними ячейками
template <int dim>
void check_link(AmrCells& cells, index_t ip, index_t main_child) {
    // Проверяем, что внутренние ячейки связаны верно
    for (int z1 = 0; z1 < CpC(dim); ++z1) {
        index_t c1 = main_child + z1;

        int count_sibs = 0;
        for (int side1: Side<dim>::items()) {
            index_t face1 = cells.face_begin[c1] + side1;

            // Грань наружу, пропускаем
            if ((cells.center[ip] - cells.center[c1]).dot(cells.faces.normal[face1]) < 0.0)
                continue;

            ++count_sibs;

            // Локальный индекс брата
            int z2 = cells.faces.adjacent.index[face1] - main_child;

            scrutiny_check(0 <= z1 && z1 < CpC(dim), "bro index in range [0, CpC(dim))")
            scrutiny_check(z1 != z2, "bro index != my index")

            // Обходим грани брата
            index_t c2 = main_child + z2;

            int side2 = 0;
            for (; side2 < Side<dim>::count(); ++side2) {
                index_t face2 = cells.face_begin[c2] + side2;

                if ((cells.faces.center[face1] - cells.faces.center[face2]).norm() < 1.0e-5 * cells.linear_size(ip)) {
                    break;
                }
            }
            if (side2 >= Side<dim>::count()) {
                std::cout << "PARENT\n";
                cells.print_info(ip);
                std::cout << "CHILD\n";
                for (int i = 0; i < CpC(dim); ++i) {
                    cells.print_info(main_child + i);
                }
            }
            // Нашли соответствующую грань, должен быть искомый
            scrutiny_check(side2 < Side<dim>::count(), "not found bro")
            index_t i3 = cells.faces.adjacent.index[cells.face_begin[c2] + side2] - main_child;
            scrutiny_check(z1 == i3, "Bad link")
        }

        // Число братьев через грань равно размерности
        if (count_sibs != dim) {
            std::cout << "Count siblings: " << count_sibs << "\n";
            std::cout << "PARENT:\n";
            cells.print_info(ip);
            std::cout << "CHILDREN:\n";
            for (int i = 0; i < CpC(dim); ++i) {
                cells.print_info(main_child + i);
            }

            check_link<dim>(cells, ip, main_child);

            throw std::runtime_error("neibs sibings count != dimension");
        }
    }
}

/// @brief Создать дочерние ячейки на выделенном месте по порядку
/// @param locals Локальное хранилище ячеек
/// @param aliens Хранилище ячеек с других процессов
/// @param ip Индекс родительской ячейки
/// @return Индекс первой дочерней ячейки, все они располагаются по порядку.
/// Дочерние ячейки имеют законченный вид (необходимое число граней, правильные
/// связи друг на друга, кроме одного случая правильные связи на соседей)
template<int dim>
index_t make_children(AmrCells &locals, AmrCells& aliens, index_t ip) {
    const index_t main_child = locals.next[ip];

    // Установить только геометрию
    if constexpr (dim == 2) {
        auto quads = locals.mapping<dim>(ip).children();
        locals.set_cell(main_child + 0, quads[0], locals.axial());
        locals.set_cell(main_child + 1, quads[1], locals.axial());
        locals.set_cell(main_child + 2, quads[2], locals.axial());
        locals.set_cell(main_child + 3, quads[3], locals.axial());
    }
    else {
        auto cubes = locals.mapping<dim>(ip).children();
        locals.set_cell(main_child + 0, cubes[0]);
        locals.set_cell(main_child + 1, cubes[1]);
        locals.set_cell(main_child + 2, cubes[2]);
        locals.set_cell(main_child + 3, cubes[3]);
        locals.set_cell(main_child + 4, cubes[4]);
        locals.set_cell(main_child + 5, cubes[5]);
        locals.set_cell(main_child + 6, cubes[6]);
        locals.set_cell(main_child + 7, cubes[7]);
    }

#if SCRUTINY
    // Бывают проблемы с выделением граней и вершин
    for (int i = 0; i < CpC(dim); ++i) {
        if (locals.max_face_count(main_child + i) != FpC(dim) * FpF(dim)) {
            throw std::runtime_error("bad max faces");
        }
        scrutiny_check(locals.max_face_count(main_child + i) == FpC(dim) * FpF(dim), "make_children error: bad faces")
        scrutiny_check(locals.max_node_count(main_child + i) == std::pow(3, dim), "make_children error: bad nodes")
    }
#endif

    int rank = locals.rank[ip];
    auto& adj = locals.faces.adjacent;

    for (int i = 0; i < CpC(dim); ++i) {
        index_t ich = main_child + i;

        locals.rank[ich] = rank;
        locals.index[ich] = ich;

        locals.flag [ich] = 0;
        locals.b_idx[ich] = locals.b_idx[ip];
        locals.level[ich] = locals.level[ip] + 1;
        locals.z_idx[ich] = CpC(dim) * locals.z_idx[ip] + i;

        // adj.basic выставим сразу у всех дочерних ячеек
        for (auto side: Side<dim>::items()) {
            adj.basic[locals.face_begin[ich] + side] = locals.next[ich];
        }
    }

    // Свяжем дочерние ячейки (выставляются boundary и все индексы смежности)
    link_siblings<dim>(locals, main_child, rank);

#if SCRUTINY
    //check_link<dim>(locals, ip, main_child);
#endif

    // Далее необходимо связать дочерние ячейки с соседями
    index_t face_beg = locals.face_begin[ip];
    for (Side<dim> side: Side<dim>::items()) {
        auto children_by_side = side.children();

        // Тип родительской грани
        auto flag = locals.faces.boundary[face_beg + side];
        if (flag == Boundary::UNDEFINED) {
            continue;
        }

        // Копируем флаг родительской на дочерние
        for (int i: children_by_side) {
            index_t ich = main_child + i;
            locals.faces.boundary[locals.face_begin[ich] + side] = flag;
        }

        // Для граничной выставляем индексы и заканчиваем
        if (locals.faces.is_boundary(face_beg + side)) {
            for (int i: children_by_side) {
                index_t ich = main_child + i;
                index_t ch_face = locals.face_begin[ich] + side;
                adj.rank [ch_face] = rank;
                adj.index[ch_face] = locals.next[ich];
                adj.alien[ch_face] = -1;
            }
            continue;
        }

        // Родительская ячейка имела простую грань по стороне
        if (locals.simple_face(ip, side)) {
            auto [neibs, jc] = adj.get_neib(face_beg + side, locals, aliens);
            scrutiny_check(0 <= jc && jc < neibs.size(), "Out fo bounds make children #1");

            scrutiny_check(
                (locals.level[ip]  > neibs.level[jc] && neibs.flag[jc] == 1) ||
                (locals.level[ip] == neibs.level[jc] && neibs.flag[jc] == 0) ||
                (locals.level[ip] == neibs.level[jc] && neibs.flag[jc] == 1),
                "All possible cases # 1");

            // Ячейка имела простую грань
            for (int i: children_by_side) {
                index_t ich = main_child + i;
                index_t ch_face = locals.face_begin[ich] + side;

                // TODO: MPI VERSION
                adj.rank[ch_face] = adj.rank[face_beg + side];
                adj.alien[ch_face] = adj.alien[face_beg + side];

                // TODO: MPI VERSION
                if (locals.level[ip] > neibs.level[jc]) {
                    // case: lvl_c > lvl_n & flag_n == 1
                    adj.index[ch_face] = locals.next[neibs.next[jc] + side.adjacent_child(locals.z_idx[ip] % CpC(dim))];
                }
                else {
                    if (neibs.flag[jc] == 0) {
                        // case: lvl_c == lvl_n & flag_n == 0
                        adj.index[ch_face] = locals.next[adj.index[face_beg + side]];
                    }
                    else {
                        // case: lvl_c == lvl_n & flag_n == 1
                        adj.index[ch_face] = locals.next[neibs.next[jc] + side.adjacent_child(i)];
                    }
                }
            }
        } else {
            // Ячейка имела сложную грань
            for (int i: children_by_side) {
                index_t ich = main_child + i;
                index_t ch_face = locals.face_begin[ich] + side;

                Side<dim> subface = side.subface_by_child(i);
                auto [neibs, jc] = adj.get_neib(face_beg + subface, locals, aliens);
                scrutiny_check(0 <= jc && jc < neibs.size(), "Out fo bounds make children #2");

                // TODO: MPI VERSION
                adj.rank[ch_face] = adj.rank[face_beg + subface];
                adj.alien[ch_face] = adj.alien[face_beg + subface];

                if (neibs.flag[jc] == 0) {
                    adj.index[ch_face] = neibs.next[jc];
                }
                else if (neibs.flag[jc] < 0) {
                    adj.index[ch_face] = locals.next[neibs.next[jc]];
                }
                else {
                    // Здесь остался возможный случай, когда сосед ещё разобьется,
                    // тогда у новой дочерней ячейки ещё придется бить грань
                    // Пока что она указывает на устаревший индекс.
                    adj.index[ch_face] = jc;
                }
            }
        }
    }

    return main_child;
}

/// @brief Производит разбиение ячейки, дочерние ячейки помещаются в хранилище
/// по порядку. Дочерние ячейки правильно ссылаются друг на друга, на гранях
/// adjacent указан правильно на новые позиции соседей.
/// @param locals Локальное хранилище ячеек
/// @param aliens Хранилище ячеек с других процессов
/// @param ip Индекс родительской ячейки
/// @param op Оператор разделения данных
template<int dim>
void refine_cell(AmrCells &locals, AmrCells& aliens, index_t ip, const Distributor& op) {
    auto main_child = make_children<dim>(locals, aliens, ip);

    auto& adj = locals.faces.adjacent;

    // Возможна ситуация, когда сторону дочерней ячейки также требуется разбить.
    // Это достигается в единственном случае, когда сосед имеет уровень выше и
    // тоже хочет разбиться. Надо проверить такие случаи.
    for (auto side: Side<dim>::items()) {
        // Интересующий нас случай встречается только тогда, когда
        // сторона родительской ячейки уже имеет подразбитую грань
        index_t p_face1 = locals.face_begin[ip] + side;
        if (locals.simple_face(ip, side) ||
            locals.faces.is_undefined(p_face1) ||
            locals.faces.is_boundary(p_face1)) {
            continue;
        }

        for (int i: side.children()) {
            index_t ich = main_child + i;
#if SCRUTINY
            if (ich < 0 || ich >= locals.size()) {
                throw std::runtime_error("[parent.next + i] out of range (refine_cell)");
            }
#endif
            index_t iface = locals.face_begin[ich] + side;
            scrutiny_check(!locals.faces.is_undefined(iface) &&
                !locals.faces.is_boundary(iface), "Should be rejected by parent check");

#if SCRUTINY
            int rank = mpi::rank();
            if (adj.rank[iface] == rank) {
                // Локальная ячейка
                if (adj.alien[iface] >= 0 || adj.index[iface] < 0 || adj.index[iface] >= locals.size()) {
                    throw std::runtime_error("adjacent.index out of range (refine_cell)");
                }
            }
            else {
                // Удаленная ячейка
                if (adj.alien[iface] < 0 || adj.alien[iface] >= aliens.size()) {
                    throw std::runtime_error("adjacent.alien out of range (refine_cell)");
                }
            }
#endif
            // Ссылка на соседнюю ячейку
            index_t p_face = locals.face_begin[ip] + side.subface_by_child(i);
            auto [neibs, jc] = adj.get_neib(p_face, locals, aliens);

            // Желаемый уровень соседней ячейки
            int neib_wanted_lvl = neibs.level[jc] + neibs.flag[jc];

            // Если сосед хочет адаптировать выше уровня текущей дочерней,
            // то дополнительно разбиваем грань дочерней ячейки
            if (neib_wanted_lvl > locals.level[ich]) {
                split_face<dim>(locals, ich, side);

                index_t neib_next = neibs.next[jc];
                for (auto subface: side.subfaces()) {
                    // TODO: MPI VERSION
                    index_t ch_face = locals.face_begin[ich] + subface;
                    adj.rank[ch_face] = adj.rank[p_face];
                    adj.index[ch_face] = locals.next[neib_next + subface.neib_child()];
                    adj.alien[ch_face] = -1;
                }
            }
        }
    }

    Children children(&locals);
    for (int i = 0; i < CpC(dim); ++i) {
        children.index[i] = main_child + i;
    }
    EuCell parent(&locals, ip);
    op.split(parent, children);

    locals.set_undefined(ip);
}

} // namespace zephyr::mesh::amr