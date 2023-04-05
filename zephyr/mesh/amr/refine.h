/// @file Файл содержит реализацию функций для разбиения ячейки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/geom/cell.h>
#include <zephyr/geom/maps.h>

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/faces.h>
#include <zephyr/mesh/amr/coarse.h>

namespace zephyr { namespace mesh { namespace amr {

/// @brief Создать геометрию дочернх ячеек по родительской ячейке
/// @param parent Родительская ячейка
/// @return Массив с дочерними ячейками
template <int dim>
std::array<geom::Cell, CpC(dim)> create_children(Storage::Item parent);

/// @brief Создать геометрию дочерних ячеек по родительской ячейке
/// @param parent Родительская ячейка
/// @return Массив с дочерними ячейками
std::array<geom::Cell, CpC(2)> create_children_simple(Storage::Item parent) {
    using geom::ShortList2D;
    using geom::Cell;

    std::array<Vector3d, 9> vs;
    vs[0] = parent.geom().vertices[0];
    vs[1] = parent.geom().vertices[1];
    vs[2] = parent.geom().vertices[2];
    vs[3] = parent.geom().vertices[3];

    vs[4] = (vs[0] + vs[2]) / 2.0;
    vs[5] = (vs[1] + vs[3]) / 2.0;
    vs[6] = (vs[0] + vs[1]) / 2.0;
    vs[7] = (vs[2] + vs[3]) / 2.0;

    vs[8] = (vs[0] + vs[1] + vs[2] + vs[3]) / 4.0;

    ShortList2D vl0 = {vs[0], vs[6], vs[4], vs[8]};
    ShortList2D vl1 = {vs[6], vs[1], vs[8], vs[5]};
    ShortList2D vl2 = {vs[4], vs[8], vs[2], vs[7]};
    ShortList2D vl3 = {vs[8], vs[5], vs[7], vs[3]};

    return {Cell(vl0), Cell(vl1), Cell(vl2), Cell(vl3)};
}

/// @brief Создать геометрию дочерних ячеек по родительским вершинам (2D)
/// @param parent_vertices Вершины родительской ячейки
/// @return Массив с дочерними ячейками
template <>
std::array<geom::Cell, CpC(2)> create_children<2>(Storage::Item parent) {
    using geom::LargeList2D;
    using geom::Mapping2D;
    using geom::Cell;

    const Vertices& vertices = parent.geom().vertices;

    // Собираем отображение ячейки
    LargeList2D vs = {
            (Vector3d &) vertices[0],
            (Vector3d &) vertices[1],
            (Vector3d &) vertices[2],
            (Vector3d &) vertices[3],
            (Vector3d &) vertices[4],
            (Vector3d &) vertices[5],
            (Vector3d &) vertices[6],
            (Vector3d &) vertices[7],
            (Vector3d &) vertices[8]
    };

    Mapping2D map(vs);
    LargeList2D vl0 = {
            vs[0], map(-0.5, -1.0), vs[1],
            map(-1.0, -0.5), map(-0.5, -0.5), map(0.0, -0.5),
            vs[3], map(-0.5, 0.0), vs[4]
    };
    LargeList2D vl1 = {
            vs[1], map(0.5, -1.0), vs[2],
            map(0.0, -0.5), map(0.5, -0.5), map(1.0, -0.5),
            vs[4], map(0.5, 0.0), vs[5]
    };
    LargeList2D vl2 = {
            vs[3], map(-0.5, 0.0), vs[4],
            map(-1.0, 0.5), map(-0.5, 0.5), map(0.0, 0.5),
            vs[6], map(-0.5, 1.0), vs[7]
    };
    LargeList2D vl3 = {
            vs[4], map(0.5, 0.0), vs[5],
            map(0.0, 0.5), map(0.5, 0.5), map(1.0, 0.5),
            vs[7], map(0.5, 1.0), vs[8]
    };

    return {Cell(vl0), Cell(vl1), Cell(vl2), Cell(vl3)};
}

/// @brief Создать геометрию дочернх ячеек по родительским вершинам (3D)
/// @param parent_vertices Вершины родительской ячейки
/// @return Массив с дочерними ячейками
template <>
std::array<Cell, CpC(3)> create_children<3>(Storage::Item parent) {
    LargeList3D& vs = (LargeList3D&) parent.geom().vertices;
   
    ShortList3D vl1 = {vs[iww(0, 0, 0)], vs[iww(1, 0, 0)], vs[iww(0, 1, 0)], vs[iww(1, 1, 0)],
                       vs[iww(0, 0, 1)], vs[iww(1, 0, 1)], vs[iww(0, 1, 1)], vs[iww(1, 1, 1)]};
    ShortList3D vl2 = {vs[iww(1, 0, 0)], vs[iww(2, 0, 0)], vs[iww(1, 1, 0)], vs[iww(2, 1, 0)],
                       vs[iww(1, 0, 1)], vs[iww(2, 0, 1)], vs[iww(1, 1, 1)], vs[iww(2, 1, 1)]};
    ShortList3D vl3 = {vs[iww(0, 1, 0)], vs[iww(1, 1, 0)], vs[iww(0, 2, 0)], vs[iww(1, 2, 0)],
                       vs[iww(0, 1, 1)], vs[iww(1, 1, 1)], vs[iww(0, 2, 1)], vs[iww(1, 2, 1)]};
    ShortList3D vl4 = {vs[iww(1, 1, 0)], vs[iww(2, 1, 0)], vs[iww(1, 2, 0)], vs[iww(2, 2, 0)],
                       vs[iww(1, 1, 1)], vs[iww(2, 1, 1)], vs[iww(1, 2, 1)], vs[iww(2, 2, 1)]};
    ShortList3D vl5 = {vs[iww(0, 0, 1)], vs[iww(1, 0, 1)], vs[iww(0, 1, 1)], vs[iww(1, 1, 1)],
                       vs[iww(0, 0, 2)], vs[iww(1, 0, 2)], vs[iww(0, 1, 2)], vs[iww(1, 1, 2)]};
    ShortList3D vl6 = {vs[iww(1, 0, 1)], vs[iww(2, 0, 1)], vs[iww(1, 1, 1)], vs[iww(2, 1, 1)],
                       vs[iww(1, 0, 2)], vs[iww(2, 0, 2)], vs[iww(1, 1, 2)], vs[iww(2, 1, 2)]};
    ShortList3D vl7 = {vs[iww(0, 1, 1)], vs[iww(1, 1, 1)], vs[iww(0, 2, 1)], vs[iww(1, 2, 1)],
                       vs[iww(0, 1, 2)], vs[iww(1, 1, 2)], vs[iww(0, 2, 2)], vs[iww(1, 2, 2)]};
    ShortList3D vl8 = {vs[iww(1, 1, 1)], vs[iww(2, 1, 1)], vs[iww(1, 2, 1)], vs[iww(2, 2, 1)],
                       vs[iww(1, 1, 2)], vs[iww(2, 1, 2)], vs[iww(1, 2, 2)], vs[iww(2, 2, 2)]};

    Vector3d zero(0.0, 0.0, 0.0);
    ShortList3D vl0 = {zero, zero, zero, zero, zero, zero, zero, zero};

    return {
            Cell(vl1), Cell(vl2), Cell(vl3), Cell(vl4),
            Cell(vl5), Cell(vl6), Cell(vl7), Cell(vl8)
    };
}

/// @brief Создать дочерние ячейки
/// @param cell Родительская ячейка
/// @param ic Индекс родительской ячейки
/// @return Массив с дочерними ячейками, дочерние ячейки имеют законченный вид
/// (необходимое число граней, правильную линковку (на старные ячейки))
template<int dim>
std::array<Cell, CpC(dim)> get_children(Storage::Item &cell, int ic, int rank) {
    auto children = create_children<dim>(cell);

    auto children_by_side = get_children_by_side<dim>();

    // По умолчанию дети ссылаются на родительскую ячейку
    for (auto &child: children) {
        for (int s = 0; s < FpC(dim); ++s) {
            child.faces[s].adjacent.rank = rank;
            child.faces[s].adjacent.index = ic;
            child.faces[s].adjacent.ghost = -1;
        }
    }

    for (int i = 0; i < CpC(dim); ++i) {
        children[i].next = cell.geom().next + i;
        children[i].flag = 0;
        children[i].b_idx = cell.geom().b_idx;
        children[i].level = cell.geom().level + 1;
        children[i].z_idx = CpC(dim) * cell.geom().z_idx + i;
    }

    // Выставить граничный флаг
    for (int side = 0; side < FpC(dim); ++side) {
        auto flag = cell.geom().faces[side].boundary;
        for (int i: children_by_side[side]) {
            children[i].faces[side].boundary = flag;
        }
    }

    // Далее необходимо слинковать дочерние ячейки с соседями
    for (int side = 0; side < FpC(dim); ++side) {
        auto flag = cell.geom().faces[side].boundary;
        for (int i: children_by_side[side]) {
            children[i].faces[side].boundary = flag;
        }

        if (cell.geom().faces[side + 6].is_undefined()) {
            // Ячейка имела простую грань
            auto adj = cell.geom().faces[side].adjacent;
            for (int i: children_by_side[side]) {
                children[i].faces[side].adjacent = adj;
            }
        } else {
            // Ячейка имела сложную грань
            for (int i: children_by_side[side]) {
                Cell &child = children[i];
                Face &child_face = child.faces[side];
                auto child_fc = child_face.center<dim>(child.vertices);

                for (auto s: subface_sides<dim>(side)) {
                    Face&cell_face = cell.geom().faces[s];
                    auto cell_fc = cell_face.center<dim>(cell.vertices());

                    if ((child_fc - cell_fc).norm() < 1.0e-5 * cell.size()) {
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
std::array<Storage::Item, CpC(dim)> select_children(
        Storage& cells, const std::array<Cell, CpC(dim)>& children);

/// @brief Возвращает итераторы дочерних ячеек (2D)
template <>
std::array<Storage::Item, CpC(2)> select_children<2>(
        Storage& cells, const std::array<Cell, CpC(2)>& children) {
    return {
            cells[children[0].next],
            cells[children[1].next],
            cells[children[2].next],
            cells[children[3].next]
    };
}

/// @brief Возвращает итераторы дочерних ячеек (3D)
template <>
std::array<Storage::Item, CpC(3)> select_children<3>(
        Storage& cells, const std::array<Cell, CpC(3)>& children) {
    return {
            cells[children[0].next],
            cells[children[1].next],
            cells[children[2].next],
            cells[children[3].next],
            cells[children[4].next],
            cells[children[5].next],
            cells[children[6].next],
            cells[children[7].next]
    };
}

/// @brief Производит разбиение ячейки, дочерние ячейки помещает в Storage
/// @param locals Хранилище ячеек
/// @param ic Индекс родительской ячейки в хранилище
/// @param op Оператор разделения данных
/// @details Дочерние ячейки правильно ссылаются друг на друга, на гранях
/// adjacent указан на старые ячейки.
template<int dim>
void refine_cell(Storage &locals, Storage &aliens, int rank, int ic, const Distributor& op) {
    auto cell = locals[ic];

    cell.set_undefined();

    auto children = get_children<dim>(cell, ic, rank);

    for (int i = 0; i < CpC(dim); ++i) {
        auto j = locals[ic].geom().next + i;

        cell.copy_to(locals[j]);

        for (int s = 0; s < FpC(dim); ++s) {
            Face &f1 = locals[j].geom().faces[s];
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
                auto neib = locals[adj.index];
                nei_wanted_lvl = neib.level() + neib.flag();
            }
            else {
                // Удаленная ячейка
#if SCRUTINY
                if (adj.ghost >= aliens.size()) {
                    throw std::runtime_error("adjacent.ghost out of range (refine_cell)");
                }
#endif
                auto neib = aliens[adj.ghost];
                nei_wanted_lvl = neib.level() + neib.flag();
            }

            if (nei_wanted_lvl > locals[j].level()) {
                split_face<dim>(locals[j], Side(s));
            }
        }
    }

    op.split<dim>(cell, select_children<dim>(locals, children));
}

} // namespace amr
} // namespace mesh
} // namespace zephyr