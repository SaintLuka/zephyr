/// @file Файл содержит реализацию функций для разбиения ячейки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/geom/maps.h>

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/faces.h>
#include <zephyr/mesh/amr/coarse.h>

namespace zephyr { namespace mesh { namespace amr {

/// @brief Создать геометрию дочернх ячеек по родительской ячейке
/// @param cube Родительская ячейка
/// @return Массив с дочерними ячейками
template <int dim>
std::array<geom::AmrCell, CpC(dim)> create_children(const SqCube& cube);

/// @brief Создать геометрию дочерних ячеек по родительским вершинам (2D)
/// @param cube Вершины родительской ячейки
/// @return Массив с дочерними ячейками
template <>
std::array<geom::AmrCell, CpC(2)> create_children<2>(const SqCube& cube) {
    auto quads = cube.as2D().children();
    return {AmrCell(quads[0]), AmrCell(quads[1]),
            AmrCell(quads[2]), AmrCell(quads[3])};
}

/// @brief Создать геометрию дочернх ячеек по родительским вершинам (3D)
/// @param cube Вершины родительской ячейки
/// @return Массив с дочерними ячейками
template <>
std::array<AmrCell, CpC(3)> create_children<3>(const SqCube& cube) {
    auto cubes = cube.children();
    return {
            AmrCell(cubes[0]), AmrCell(cubes[1]), AmrCell(cubes[2]), AmrCell(cubes[3]),
            AmrCell(cubes[4]), AmrCell(cubes[5]), AmrCell(cubes[6]), AmrCell(cubes[7])
    };
}

/// @brief Создать дочерние ячейки
/// @param cell Родительская ячейка
/// @param ic Индекс родительской ячейки
/// @return Массив с дочерними ячейками, дочерние ячейки имеют законченный вид
/// (необходимое число граней, правильную линковку (на старные ячейки))
template<int dim>
std::array<AmrCell, CpC(dim)> get_children(AmrCell &cell) {
    auto children = create_children<dim>(cell.vertices);

    auto children_by_side = get_children_by_side<dim>();

    // По умолчанию дети ссылаются на родительскую ячейку
    for (auto &child: children) {
        child.rank = cell.rank;
        child.index = cell.index;

        for (int s = 0; s < FpC(dim); ++s) {
            child.faces[s].adjacent.rank = cell.rank;
            child.faces[s].adjacent.index = cell.index;
            child.faces[s].adjacent.ghost = -1;
        }
    }

    for (int i = 0; i < CpC(dim); ++i) {
        children[i].next = cell.next + i;
        children[i].flag = 0;
        children[i].b_idx = cell.b_idx;
        children[i].level = cell.level + 1;
        children[i].z_idx = CpC(dim) * cell.z_idx + i;
    }

    // Выставить граничный флаг
    for (int side = 0; side < FpC(dim); ++side) {
        auto flag = cell.faces[side].boundary;
        for (int i: children_by_side[side]) {
            children[i].faces[side].boundary = flag;
        }
    }

    // Далее необходимо слинковать дочерние ячейки с соседями
    for (int side = 0; side < FpC(dim); ++side) {
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
                AmrCell &child = children[i];
                BFace &child_face = child.faces[side];
                auto child_fc = child_face.center;

                for (auto s: subface_sides<dim>(side)) {
                    BFace& cell_face = cell.faces[s];
                    auto cell_fc = cell_face.center;

                    if ((child_fc - cell_fc).norm() < 1.0e-5 * cell.size) {
                        child_face.adjacent = cell_face.adjacent;
                        break;
                    }
                }
            }
        }
    }

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
/// @param locals Хранилище ячеек
/// @param ic Индекс родительской ячейки в хранилище
/// @param op Оператор разделения данных
/// @details Дочерние ячейки правильно ссылаются друг на друга, на гранях
/// adjacent указан на старые ячейки.
template<int dim>
void refine_cell(AmrStorage &locals, AmrStorage &aliens, int rank, int ic, const Distributor& op) {
    auto& parent = locals[ic];

    auto children = get_children<dim>(parent);

    for (int i = 0; i < CpC(dim); ++i) {
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
                if (adj.index >= locals.size()) {
                    throw std::runtime_error("adjacent.index out of range (refine_cell)");
                }
#endif
                auto& neib = locals[adj.index];
                nei_wanted_lvl = neib.level + neib.flag;
            }
            else {
                // Удаленная ячейка
#if SCRUTINY
                if (adj.ghost >= aliens.size()) {
                    throw std::runtime_error("adjacent.ghost out of range (refine_cell)");
                }
#endif
                auto& neib = aliens[adj.ghost];
                nei_wanted_lvl = neib.level + neib.flag;
            }

            if (nei_wanted_lvl > child.level) {
                split_face<dim>(child, Side(s));
            }
        }
    }

    auto children2 = select_children<dim>(locals, children);
    op.split(parent, children2);

    parent.set_undefined();
}

} // namespace amr
} // namespace mesh
} // namespace zephyr