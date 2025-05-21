/// @brief Файл содержит реализацию функций для разбиения ячейки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/primitives/amr_cell.h>
#include <zephyr/geom/primitives/cube.h>

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/faces.h>
#include <zephyr/mesh/amr/coarse.h>

namespace zephyr::mesh::amr {

/// @brief Создать геометрию дочернх ячеек по родительской ячейке
/// @param cube Родительская ячейка
/// @return Массив с дочерними ячейками
template <int dim>
std::array<AmrCell, CpC(dim)> create_children(const SqCube& cube, bool axial);

/// @brief Создать геометрию дочерних ячеек по родительским вершинам (2D)
/// @param cube Вершины родительской ячейки
/// @return Массив с дочерними ячейками
template <>
std::array<AmrCell, CpC(2)> create_children<2>(const SqCube& cube, bool axial) {
    auto quads = cube.as2D().children();
    return {AmrCell(quads[0], axial), AmrCell(quads[1], axial),
            AmrCell(quads[2], axial), AmrCell(quads[3], axial)};
}

/// @brief Создать геометрию дочернх ячеек по родительским вершинам (3D)
/// @param cube Вершины родительской ячейки
/// @return Массив с дочерними ячейками
template <>
std::array<AmrCell, CpC(3)> create_children<3>(const SqCube& cube, bool axial) {
    auto cubes = cube.children();
    return {
            AmrCell(cubes[0]), AmrCell(cubes[1]), AmrCell(cubes[2]), AmrCell(cubes[3]),
            AmrCell(cubes[4]), AmrCell(cubes[5]), AmrCell(cubes[6]), AmrCell(cubes[7])
    };
}

/// @brief Связать грани соседних дочерних ячеек
/// @param children Набор дочерних ячеек
/// @param first_index Индекс первой ячейки
/// Помним, что все дочерние ячейки будут располагаться последовательно,
/// нумерация начнется с first_index.
template <int dim>
void link_siblings(std::array<AmrCell, CpC(dim)>& children, int first_index);

#define subs2D(i, j, x, y) (children[Quad::iss<i, j>()].faces[Side2D::by_dir<x, y>()].adjacent.index = first_index + Quad::iss<i + 2 * x, j + 2 * y>())
#define subs3D(i, j, k, x, y, z) (children[Cube::iss<i, j, k>()].faces[Side3D::by_dir<x, y, z>()].adjacent.index = first_index + Cube::iss<i + 2 * x, j + 2 * y, k + 2 * z>())

template <>
inline void link_siblings<2>(std::array<AmrCell, CpC(2)>& children, int first_index) {
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
inline void link_siblings<3>(std::array<AmrCell, CpC(3)>& children, int first_index) {
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
void check_link(std::array<AmrCell, CpC(dim)>& children, AmrCell &parent) {
    // Проверяем, что внутренние ячейки связаны верно
    for (int c1 = 0; c1 < CpC(dim); ++c1) {
        auto& child1 = children[c1];

        int count_sibs = 0;
        for (int side1 = 0; side1 < FpC(dim); ++side1) {
            auto& face1 = child1.faces[side1];

            // Грань наружу, пропускаем
            if ((parent.center - child1.center).dot(face1.normal) < 0.0)
                continue;

            ++count_sibs;

            // Локальный индекс брата
            int c2 = face1.adjacent.index - parent.next;

            scrutiny_check(0 <= c1 && c1 < CpC(dim), "bro index in range [0, CpC(dim))")
            scrutiny_check(c1 != c2, "bro index != my index")

            // Обходим грани брата
            auto& child2 = children[c2];

            int side2 = 0;
            for (; side2 < FpC(dim); ++side2) {
                auto &face2 = child2.faces[side2];

                if ((face1.center - face2.center).norm() < 1.0e-5 * parent.linear_size()) {
                    break;
                }
            }
            // Нашли соответствующую грань, должен быть искомый
            scrutiny_check(side2 < FpC(dim), "not found bro")
            int c3 = child2.faces[side2].adjacent.index - parent.next;
            scrutiny_check(c1 == c3, "Bad link")
        }

        // Число братьев через грань равно размерности
        if (count_sibs != dim) {
            std::cout << "Count siblings: " << count_sibs << "\n";
            throw std::runtime_error("neibs sibings count != dimension");
        }
    }
}

/// @brief Создать дочерние ячейки
/// @param cell Родительская ячейка
/// @param ic Индекс родительской ячейки
/// @return Массив с дочерними ячейками, дочерние ячейки имеют законченный вид
/// (необходимое число граней, правильную линковку (на старные ячейки))
template<int dim>
std::array<AmrCell, CpC(dim)> get_children(AmrCell &cell) {
    const auto children_by_side = get_children_by_side<dim>();

    auto children = create_children<dim>(cell.vertices, cell.axial);

    for (int i = 0; i < CpC(dim); ++i) {
        auto& child = children[i];

        child.rank  = cell.rank;
        child.index = cell.next + i;

        child.next = cell.next + i;
        child.flag = 0;
        child.b_idx = cell.b_idx;
        child.level = cell.level + 1;
        child.z_idx = CpC(dim) * cell.z_idx + i;

        // По умолчанию дети ссылаются на родительскую ячейку
        // зачем этот код? Потом все изменяются
        for (int s = 0; s < FpC(dim); ++s) {
            child.faces[s].adjacent.rank  = cell.rank;
            child.faces[s].adjacent.index = cell.index;
            child.faces[s].adjacent.alien = -1;
        }
    }

    // Далее необходимо связать дочерние ячейки с соседями
    for (int side = 0; side < FpC(dim); ++side) {
        // Выставить граничный флаг
        auto flag = cell.faces[side].boundary;
        for (int i: children_by_side[side]) {
            children[i].faces[side].boundary = flag;
        }

        if (cell.faces[side + 6].is_undefined()) {
            // Ячейка имела простую грань
            auto adj = cell.faces[side].adjacent;
            for (int i: children_by_side[side]) {
                children[i].faces[side].adjacent = adj;
            }
        } else {
            // Ячейка имела сложную грань
            for (int i: children_by_side[side]) {
                BFace &child_face = children[i].faces[side];
                auto child_fc = child_face.center;

                for (auto s: subface_sides<dim>(side)) {
                    BFace& cell_face = cell.faces[s];
                    auto cell_fc = cell_face.center;

                    if ((child_fc - cell_fc).norm() < 1.0e-5 * cell.linear_size()) {
                        child_face.adjacent = cell_face.adjacent;
                        break;
                    }
                }
            }
        }
    }

    // Свяжем внутренние ячейки
    link_siblings<dim>(children, cell.next);

#if SCRUTINY
    check_link<dim>(children, cell);
#endif

    return children;
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
/// @param locals Хранилище ячеек
/// @param op Оператор разделения данных
/// @details Дочерние ячейки правильно ссылаются друг на друга, на гранях
/// adjacent указан на старые ячейки.
template<int dim>
void refine_cell(AmrStorage::Item& parent, AmrStorage &locals, AmrStorage &aliens, int rank, const Distributor& op) {
    auto children = get_children<dim>(parent);

    for (int i = 0; i < CpC(dim); ++i) {
#if SCRUTINY
        if (parent.next + i < 0 || parent.next + i >= locals.size()) {
            throw std::runtime_error("[parent.next + i] out of range (refine_cell)");
        }
#endif

        AmrCell& child = locals[parent.next + i];

        // Возможна оптимизация (!), создавать дочерние ячейки сразу
        // на нужном месте в хранилище
        locals[parent.next + i] = children[i];

        for (int s = 0; s < FpC(dim); ++s) {
            BFace &f1 = child.faces[s];
            if (f1.is_undefined() or f1.is_boundary()) {
                continue;
            }

            auto adj = f1.adjacent;
            int nei_wanted_lvl = 0;
            if (adj.rank == rank) {
                // Локальная ячейка
#if SCRUTINY
                if (adj.index < 0 || adj.index >= locals.size()) {
                    throw std::runtime_error("adjacent.index out of range (refine_cell)");
                }
#endif
                auto& neib = locals[adj.index];
                nei_wanted_lvl = neib.level + neib.flag;
            }
            else {
                // Удаленная ячейка
#if SCRUTINY
                if (adj.alien < 0 || adj.alien >= aliens.size()) {
                    throw std::runtime_error("adjacent.alien out of range (refine_cell)");
                }
#endif
                auto& neib = aliens[adj.alien];
                nei_wanted_lvl = neib.level + neib.flag;
            }

            if (nei_wanted_lvl > child.level) {
                split_face<dim>(child, Side3D(s));
            }
        }
    }

    auto children2 = select_children<dim>(locals, children);
    op.split(parent, children2);

    parent.set_undefined();
}

} // namespace zephyr::mesh::amr