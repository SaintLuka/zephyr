/// @brief Файл содержит реализацию функций для разбиения ячейки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/primitives/amr_cell.h>
#include <zephyr/geom/primitives/cube.h>

#include <zephyr/mesh/amr2/common.h>
#include <zephyr/mesh/amr2/faces.h>
#include <zephyr/mesh/amr2/coarse.h>

namespace zephyr::mesh::amr2 {

/// @brief Создать геометрию дочернх ячеек по родительской ячейке
/// @param ic Позиция для размещения дочерних
/// @param shape Родительская ячейка
template <int dim>
void create_children(index_t ic, SoaCell& cells, const SqMap<dim>& shape);

/// @brief Создать геометрию дочерних ячеек по родительским вершинам (2D)
/// @param quad Вершины родительской ячейки
template <>
inline void create_children<2>(index_t ic, SoaCell& cells, const SqQuad& quad) {
    auto quads = quad.children();
    cells.add_cell(ic + 0, quads[0], cells.axial);
    cells.add_cell(ic + 1, quads[1], cells.axial);
    cells.add_cell(ic + 2, quads[2], cells.axial);
    cells.add_cell(ic + 3, quads[3], cells.axial);
}

/// @brief Создать геометрию дочернх ячеек по родительским вершинам (3D)
/// @param cube Вершины родительской ячейки
template <>
inline void create_children<3>(index_t ic, SoaCell& cells, const SqCube& cube) {
    auto cubes = cube.children();
    cells.add_cell(ic + 0, cubes[0]);
    cells.add_cell(ic + 1, cubes[1]);
    cells.add_cell(ic + 2, cubes[2]);
    cells.add_cell(ic + 3, cubes[3]);
    cells.add_cell(ic + 4, cubes[4]);
    cells.add_cell(ic + 5, cubes[5]);
    cells.add_cell(ic + 6, cubes[6]);
    cells.add_cell(ic + 7, cubes[7]);
}

/// @brief Связать грани соседних дочерних ячеек
/// @param ic Индекс первой ячейки
/// Помним, что все дочерние ячейки будут располагаться последовательно,
/// нумерация начнется с ic.
template <int dim>
void link_siblings(index_t ic, SoaCell& cells);

// OLD VERSION
//#define subs2D(i, j, x, y) (children[Quad::iss<i, j>()].faces[side_by_dir<x, y>()].adjacent.index = first_index + Quad::iss<i + 2 * x, j + 2 * y>())
//#define subs3D(i, j, k, x, y, z) (children[Cube::iss<i, j, k>()].faces[side_by_dir<x, y, z>()].adjacent.index = first_index + Cube::iss<i + 2 * x, j + 2 * y, k + 2 * z>())

#define subs2D(i, j, x, y) (cells.faces.adjacent.local_index[cells.face_begin[ic] + Quad::iss<i, j>()] = ic + Quad::iss<i + 2 * x, j + 2 * y>())
#define subs3D(i, j, k, x, y, z) (cells.faces.adjacent.local_index[cells.face_begin[ic] + Cube::iss<i, j, k>()] = ic + Cube::iss<i + 2 * x, j + 2 * y, k + 2 * z>())

template <>
inline void link_siblings<2>(index_t ic, SoaCell& cells) {
    subs2D(-1, -1, +1, 0);
    subs2D(-1, -1, 0, +1);

    subs2D(+1, -1, -1, 0);
    subs2D(+1, -1, 0, +1);

    subs2D(-1, +1, +1, 0);
    subs2D(-1, +1, 0, -1);

    subs2D(+1, +1, -1, 0);
    subs2D(+1, +1, 0, -1);
}

template <>
inline void link_siblings<3>(index_t ic, SoaCell& cells) {
    subs3D(-1, -1, -1, +1, 0, 0);
    subs3D(-1, -1, -1, 0, +1, 0);
    subs3D(-1, -1, -1, 0, 0, +1);

    subs3D(+1, -1, -1, -1, 0, 0);
    subs3D(+1, -1, -1, 0, +1, 0);
    subs3D(+1, -1, -1, 0, 0, +1);

    subs3D(-1, +1, -1, +1, 0, 0);
    subs3D(-1, +1, -1, 0, -1, 0);
    subs3D(-1, +1, -1, 0, 0, +1);

    subs3D(+1, +1, -1, -1, 0, 0);
    subs3D(+1, +1, -1, 0, -1, 0);
    subs3D(+1, +1, -1, 0, 0, +1);

    subs3D(-1, -1, +1, +1, 0, 0);
    subs3D(-1, -1, +1, 0, +1, 0);
    subs3D(-1, -1, +1, 0, 0, -1);

    subs3D(+1, -1, +1, -1, 0, 0);
    subs3D(+1, -1, +1, 0, +1, 0);
    subs3D(+1, -1, +1, 0, 0, -1);

    subs3D(-1, +1, +1, +1, 0, 0);
    subs3D(-1, +1, +1, 0, -1, 0);
    subs3D(-1, +1, +1, 0, 0, -1);

    subs3D(+1, +1, +1, -1, 0, 0);
    subs3D(+1, +1, +1, 0, -1, 0);
    subs3D(+1, +1, +1, 0, 0, -1);
}

template <int dim>
void check_link(index_t ip, index_t main_child, SoaCell& cells) {
    // Проверяем, что внутренние ячейки связаны верно
    for (int i1 = 0; i1 < CpC(dim); ++i1) {
        index_t c1 = main_child + i1;

        int count_sibs = 0;
        for (int side1 = 0; side1 < FpC(dim); ++side1) {
            index_t face1 = cells.face_begin[c1] + side1;

            // Грань наружу, пропускаем
            if ((cells.center[ip] - cells.center[c1]).dot(cells.faces.normal[face1]) < 0.0)
                continue;

            ++count_sibs;

            // Локальный индекс брата
            int i2 = cells.faces.adjacent.local_index[face1] - main_child;

            scrutiny_check(0 <= i1 && i1 < CpC(dim), "bro index in range [0, CpC(dim))")
            scrutiny_check(i1 != i2, "bro index != my index")

            // Обходим грани брата
            index_t c2 = main_child + i2;

            int side2 = 0;
            for (; side2 < FpC(dim); ++side2) {
                index_t face2 = cells.face_begin[c2] + side2;

                if ((cells.faces.center[face1] - cells.faces.center[face2]).norm() < 1.0e-5 * cells.linear_size(ip)) {
                    break;
                }
            }
            // Нашли соответствующую грань, должен быть искомый
            scrutiny_check(side2 < FpC(dim), "not found bro")
            index_t i3 = cells.faces.adjacent.local_index[cells.face_begin[c2] + side2] - main_child;
            scrutiny_check(i1 == i3, "Bad link")
        }

        // Число братьев через грань равно размерности
        if (count_sibs != dim) {
            std::cout << "Count siblings: " << count_sibs << "\n";
            throw std::runtime_error("neibs sibings count != dimension");
        }
    }
}

/// @brief Создать дочерние ячейки на выделенном месте по порядку
/// @return Массив индексов дочерних ячеек, дочерние ячейки имеют законченный вид
/// (необходимое число граней, правильную линковку (на старные ячейки))
template<int dim>
index_t make_children(index_t ip, SoaCell &cells) {
    const auto children_by_side = get_children_by_side<dim>();

    auto& adj = cells.faces.adjacent;

    index_t main_child = cells.next[ip];
    index_t p_face = cells.face_begin[ip];

    create_children<dim>(main_child, cells, cells.get_vertices<dim>(ip));

    for (int i = 0; i < CpC(dim); ++i) {
        index_t ich = main_child + i;

        cells.rank[ich]  = cells.rank[ip];
        cells.owner_index[ich] = ich;

        cells.next[ich] = ich;
        cells.flag[ich] = 0;
        cells.b_idx[ich] = cells.b_idx[ip];
        cells.level[ich] = cells.level[ip] + 1;
        cells.z_idx[ich] = CpC(dim) * cells.z_idx[ip] + i;

        // По умолчанию дети ссылаются на родительскую ячейку
        // зачем этот код? Потом все изменяются
        for (int s = 0; s < FpC(dim); ++s) {
            adj.rank[cells.face_begin[ich] + s]        = cells.rank[ip];
            adj.owner_index[cells.face_begin[ich] + s] = cells.owner_index[ip];
            adj.local_index[cells.face_begin[ich] + s] = cells.owner_index[ip];
        }
    }

    // Далее необходимо связать дочерние ячейки с соседями
    for (int side = 0; side < FpC(dim); ++side) {
        // Выставить граничный флаг
        auto flag = cells.faces.boundary[p_face + side];
        for (int i: children_by_side[side]) {
            index_t ich = main_child + i;
            cells.faces.boundary[cells.face_begin[ich] + side] = flag;
        }

        if (cells.faces.is_undefined(p_face + side + 6)) {
            // Ячейка имела простую грань
            for (int i: children_by_side[side]) {
                index_t ich = main_child + i;
                index_t ch_face = cells.face_begin[ich] + side;

                adj.rank[ch_face]        = adj.rank[p_face + side];
                adj.owner_index[ch_face] = adj.owner_index[p_face + side];
                adj.local_index[ch_face] = adj.local_index[p_face + side];
            }
        } else {
            // Ячейка имела сложную грань
            for (int i: children_by_side[side]) {
                index_t ich = main_child + i;
                index_t ch_face = cells.face_begin[ich] + side;

                auto child_fc = cells.faces.center[ch_face];

                for (auto s: subface_sides<dim>(side)) {
                    auto cell_fc = cells.faces.center[p_face];

                    if ((child_fc - cell_fc).norm() < 1.0e-5 * cells.linear_size(ip)) {
                        adj.rank[ch_face]        = adj.rank[p_face + side];
                        adj.owner_index[ch_face] = adj.owner_index[p_face + side];
                        adj.local_index[ch_face] = adj.local_index[p_face + side];
                        break;
                    }
                }
            }
        }
    }

    // Свяжем внутренние ячейки
    link_siblings<dim>(main_child, cells);

#if SCRUTINY
    check_link<dim>(ip, main_child, cells);
#endif

    return main_child;
}

/// @brief Возвращает итераторы дочерних ячеек
template <int dim>
Children select_children(
        AmrStorage& cells, const std::array<AmrCell, CpC(dim)>& children);

/// @brief Возвращает итераторы дочерних ячеек (2D)
template <>
Children select_children<2>(
        AmrStorage& cells, const std::array<AmrCell, CpC(2)>& children) {
    return {
            cells.iterator(children[0].next),
            cells.iterator(children[1].next),
            cells.iterator(children[2].next),
            cells.iterator(children[3].next)
    };
}

/// @brief Возвращает итераторы дочерних ячеек (3D)
template <>
Children select_children<3>(
        AmrStorage& cells, const std::array<AmrCell, CpC(3)>& children) {
    return {
            cells.iterator(children[0].next),
            cells.iterator(children[1].next),
            cells.iterator(children[2].next),
            cells.iterator(children[3].next),
            cells.iterator(children[4].next),
            cells.iterator(children[5].next),
            cells.iterator(children[6].next),
            cells.iterator(children[7].next)
    };
}

/// @brief Производит разбиение ячейки, дочерние ячейки помещает в AmrStorage
/// @param parent Родительская ячейка в хранилище
/// @param cells Хранилище ячеек
/// @param op Оператор разделения данных
/// @details Дочерние ячейки правильно ссылаются друг на друга, на гранях
/// adjacent указан на старые ячейки.
template<int dim>
void refine_cell(index_t ip, SoaCell &cells, int rank, const Distributor& op) {
    auto main_child = make_children<dim>(ip, cells);

    auto& adj = cells.faces.adjacent;

    for (int i = 0; i < CpC(dim); ++i) {
#if SCRUTINY
        if (main_child + i < 0 || main_child + i >= cells.n_cells()) {
            throw std::runtime_error("[parent.next + i] out of range (refine_cell)");
        }
#endif
        index_t ich = main_child + i;

        for (int s = 0; s < FpC(dim); ++s) {
            index_t iface = cells.face_begin[ich] + s;
            if (cells.faces.is_undefined(iface) or
                cells.faces.is_boundary(iface)) {
                continue;
            }

            int nei_wanted_lvl = 0;
            if (adj.rank[iface] == rank) {
                // Локальная ячейка
#if SCRUTINY
                if (adj.local_index[iface] < 0 || adj.local_index[iface] >= cells.n_cells()) {
                    throw std::runtime_error("adjacent.index out of range (refine_cell)");
                }
#endif
                index_t neib_idx = adj.local_index[iface];
                nei_wanted_lvl = cells.level[neib_idx] + cells.flag[neib_idx];
            }
            else {
                // Удаленная ячейка
#if SCRUTINY
                if (adj.owner_index[iface] < 0 || adj.local_index[iface] >= cells.n_cells()) {
                    throw std::runtime_error("adjacent.alien out of range (refine_cell)");
                }
#endif
                index_t neib_idx = adj.local_index[iface];
                nei_wanted_lvl = cells.level[neib_idx] + cells.flag[neib_idx];
            }

            if (nei_wanted_lvl > cells.level[ich]) {
                split_face<dim>(iface, cells.faces, cells.get_vertices<dim>(ich), s, cells.axial);
            }
        }
    }

    //auto children2 = select_children<dim>(cells, children);
    //op.split(parent, children2);

    cells.set_undefined(ip);
}

} // namespace zephyr::mesh::amr2