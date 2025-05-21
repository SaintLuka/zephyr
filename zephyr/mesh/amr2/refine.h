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

/// @brief Создать геометрию дочерних ячеек по родительской ячейке
/// @param ic Позиция для размещения дочерних
/// @param shape Родительская ячейка (SqQuad или SqCube
template <int dim>
void create_children(index_t ic, SoaCell& cells, const SqMap<dim>& shape) {
    if constexpr (dim == 2) {
        auto quads = shape.children();
        cells.add_cell(ic + 0, quads[0], cells.axial);
        cells.add_cell(ic + 1, quads[1], cells.axial);
        cells.add_cell(ic + 2, quads[2], cells.axial);
        cells.add_cell(ic + 3, quads[3], cells.axial);
    }
    else {
        auto cubes = shape.children();
        cells.add_cell(ic + 0, cubes[0]);
        cells.add_cell(ic + 1, cubes[1]);
        cells.add_cell(ic + 2, cubes[2]);
        cells.add_cell(ic + 3, cubes[3]);
        cells.add_cell(ic + 4, cubes[4]);
        cells.add_cell(ic + 5, cubes[5]);
        cells.add_cell(ic + 6, cubes[6]);
        cells.add_cell(ic + 7, cubes[7]);
    }
}

/// @brief Связать грани соседних дочерних ячеек
/// @param ic Индекс первой ячейки
/// Помним, что все дочерние ячейки будут располагаться последовательно,
/// нумерация начнется с ic.
template <int dim>
void link_siblings(SoaCell& cells, index_t ic);

#define subs2D(i, j, x, y) (cells.faces.adjacent.index[cells.face_begin[ic + Quad::iss<i, j>()] + Side2D::by_dir<x, y>()] = ic + Quad::iss<i + 2 * x, j + 2 * y>())
#define subs3D(i, j, k, x, y, z) (cells.faces.adjacent.index[cells.face_begin[ic + Cube::iss<i, j, k>()] +  + Side3D::by_dir<x, y, z>()] = ic + Cube::iss<i + 2 * x, j + 2 * y, k + 2 * z>())

template <>
inline void link_siblings<2>(SoaCell& cells, index_t ic) {
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
inline void link_siblings<3>(SoaCell& cells, index_t ic) {
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
void check_link(SoaCell& cells, index_t ip, index_t child_beg) {
    // Проверяем, что внутренние ячейки связаны верно
    for (int z1 = 0; z1 < CpC(dim); ++z1) {
        index_t c1 = child_beg + z1;

        int count_sibs = 0;
        for (int side1: Side<dim>::items()) {
            index_t face1 = cells.face_begin[c1] + side1;

            // Грань наружу, пропускаем
            if ((cells.center[ip] - cells.center[c1]).dot(cells.faces.normal[face1]) < 0.0)
                continue;

            ++count_sibs;

            // Локальный индекс брата
            int z2 = cells.faces.adjacent.index[face1] - child_beg;

            scrutiny_check(0 <= z1 && z1 < CpC(dim), "bro index in range [0, CpC(dim))")
            scrutiny_check(z1 != z2, "bro index != my index")

            // Обходим грани брата
            index_t c2 = child_beg + z2;

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
                    cells.print_info(child_beg + i);
                }
            }
            // Нашли соответствующую грань, должен быть искомый
            scrutiny_check(side2 < Side<dim>::count(), "not found bro")
            index_t i3 = cells.faces.adjacent.index[cells.face_begin[c2] + side2] - child_beg;
            scrutiny_check(z1 == i3, "Bad link")
        }

        // Число братьев через грань равно размерности
        if (count_sibs != dim) {
            std::cout << "Count siblings: " << count_sibs << "\n";
            throw std::runtime_error("neibs sibings count != dimension");
        }
    }
}

/// @brief Создать дочерние ячейки на выделенном месте по порядку
/// @param cells Хранилище ячеек
/// @param ip Индекс родительской ячейки
/// @return Индекс первой дочерней ячейки, все они располагаются по порядку.
/// Дочерние ячейки имеют законченный вид (необходимое число граней,
/// правильную линковку на старые ячейки)
template<int dim>
index_t make_children(SoaCell &cells, index_t ip) {
    const auto children_by_side = get_children_by_side<dim>();

    const index_t child_beg = cells.next[ip];

    if constexpr (dim == 2) {
        auto quads = cells.get_vertices<dim>(ip).children();
        cells.add_cell(child_beg + 0, quads[0], cells.axial);
        cells.add_cell(child_beg + 1, quads[1], cells.axial);
        cells.add_cell(child_beg + 2, quads[2], cells.axial);
        cells.add_cell(child_beg + 3, quads[3], cells.axial);
    }
    else {
        auto cubes = cells.get_vertices<dim>(ip).children();
        cells.add_cell(child_beg + 0, cubes[0]);
        cells.add_cell(child_beg + 1, cubes[1]);
        cells.add_cell(child_beg + 2, cubes[2]);
        cells.add_cell(child_beg + 3, cubes[3]);
        cells.add_cell(child_beg + 4, cubes[4]);
        cells.add_cell(child_beg + 5, cubes[5]);
        cells.add_cell(child_beg + 6, cubes[6]);
        cells.add_cell(child_beg + 7, cubes[7]);
    }

#ifdef SCRUTINY
    // Бывают проблемы с выделением граней и вершин
    for (int i = 0; i < CpC(dim); ++i) {
        if (cells.max_faces(child_beg + i) != FpC(dim) * FpF(dim)) {

            throw std::runtime_error("bad max faces");
        }
        scrutiny_check(cells.max_faces(child_beg + i) == FpC(dim) * FpF(dim), "make_children error: bad faces")
        scrutiny_check(cells.max_nodes(child_beg + i) == std::pow(3, dim), "make_children error: bad nodes")
    }
#endif

    auto& adj = cells.faces.adjacent;

    index_t p_face = cells.face_begin[ip];

    for (int i = 0; i < CpC(dim); ++i) {
        index_t ich = child_beg + i;

        cells.rank[ich]  = cells.rank[ip];
        cells.index[ich] = ich;

        cells.next[ich] = ich;
        cells.flag[ich] = 0;
        cells.b_idx[ich] = cells.b_idx[ip];
        cells.level[ich] = cells.level[ip] + 1;
        cells.z_idx[ich] = CpC(dim) * cells.z_idx[ip] + i;

        // По умолчанию дети ссылаются на родительскую ячейку
        // зачем этот код? Потом все изменяются
        for (auto s: Side<dim>::items()) {
            adj.rank[cells.face_begin[ich] + s]  = cells.rank[ip];
            adj.index[cells.face_begin[ich] + s] = cells.index[ip];
            adj.alien[cells.face_begin[ich] + s] = -1;
        }
    }

    // Далее необходимо связать дочерние ячейки с соседями
    for (Side<dim> side: Side<dim>::items()) {
        // Выставить граничный флаг
        auto flag = cells.faces.boundary[p_face + side];
        for (int i: children_by_side[side]) {
            index_t ich = child_beg + i;
            cells.faces.boundary[cells.face_begin[ich] + side] = flag;
        }

        if (cells.simple_face(ip, side)) {
            // Ячейка имела простую грань
            for (int i: children_by_side[side]) {
                index_t ich = child_beg + i;
                index_t ch_face = cells.face_begin[ich] + side;

                adj.rank [ch_face] = adj.rank [p_face + side];
                adj.index[ch_face] = adj.index[p_face + side];
                adj.alien[ch_face] = adj.alien[p_face + side];
            }
        } else {
            // Ячейка имела сложную грань
            for (int i: children_by_side[side]) {
                index_t ich = child_beg + i;
                index_t ch_face = cells.face_begin[ich] + side;

                auto child_fc = cells.faces.center[ch_face];

                for (auto s: subface_sides<dim>(side)) {
                    auto cell_fc = cells.faces.center[p_face + s];

                    if ((child_fc - cell_fc).norm() < 1.0e-5 * cells.linear_size(ip)) {
                        adj.rank [ch_face] = adj.rank [p_face + s];
                        adj.index[ch_face] = adj.index[p_face + s];
                        adj.alien[ch_face] = adj.alien[p_face + s];
                        break;
                    }
                }
            }
        }
    }

    // Свяжем внутренние ячейки
    link_siblings<dim>(cells, child_beg);

#if SCRUTINY
    check_link<dim>(cells, ip, child_beg);
#endif

    return child_beg;
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

/// @brief Производит разбиение ячейки, дочерние ячейки помещает в хранилище
/// @param locals Хранилище ячеек
/// @param ip Индекс родительской ячейки
/// @param op Оператор разделения данных
/// @details Дочерние ячейки правильно ссылаются друг на друга, на гранях
/// adjacent указан на старые ячейки.
template<int dim>
void refine_cell(SoaCell &locals, index_t ip, int rank, const Distributor& op) {
    auto child_beg = make_children<dim>(locals, ip);

    auto& adj = locals.faces.adjacent;

    for (int i = 0; i < CpC(dim); ++i) {
#if SCRUTINY
        if (child_beg + i < 0 || child_beg + i >= locals.size()) {
            throw std::runtime_error("[parent.next + i] out of range (refine_cell)");
        }
#endif
        index_t ich = child_beg + i;

        for (Side<dim> s: Side<dim>::items()) {
            index_t iface = locals.face_begin[ich] + s;
            if (locals.faces.is_undefined(iface) or
                locals.faces.is_boundary(iface)) {
                continue;
            }

            int nei_wanted_lvl = 0;
            if (adj.rank[iface] == rank) {
                // Локальная ячейка
#if SCRUTINY
                if (adj.index[iface] < 0 || adj.index[iface] >= locals.size()) {
                    throw std::runtime_error("adjacent.index out of range (refine_cell)");
                }
#endif
                index_t neib_idx = adj.index[iface];
                nei_wanted_lvl = locals.level[neib_idx] + locals.flag[neib_idx];
            }
            else {
                // Удаленная ячейка
#if SCRUTINY
                if (adj.alien[iface] < 0 || adj.index[iface] >= locals.size()) {
                    throw std::runtime_error("adjacent.alien out of range (refine_cell)");
                }
#endif
                index_t neib_idx = adj.index[iface];
                nei_wanted_lvl = locals.level[neib_idx] + locals.flag[neib_idx];
            }

            if (nei_wanted_lvl > locals.level[ich]) {
                split_face<dim>(locals.face_begin[ich], locals.faces,
                        locals.get_vertices<dim>(ich), s, locals.axial);
            }
        }
    }

    //auto children2 = select_children<dim>(cells, children);
    //op.split(parent, children2);

    locals.set_undefined(ip);
}

} // namespace zephyr::mesh::amr2