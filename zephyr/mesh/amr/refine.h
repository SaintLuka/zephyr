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
/// @param parent Родительская ячейка
/// @return Массив с дочерними ячейками
template <int dim>
std::array<geom::AmrCell, CpC(dim)> create_children(AmrCell& parent);

/// @brief Создать геометрию дочерних ячеек по родительской ячейке
/// @param parent Родительская ячейка
/// @return Массив с дочерними ячейками
std::array<geom::AmrCell, CpC(2)> create_children_simple(AmrCell& parent) {
    using geom::Quad;
    using geom::AmrCell;

    std::array<Vector3d, 9> vs;
    vs[0] = parent.vertices[0];
    vs[1] = parent.vertices[1];
    vs[2] = parent.vertices[2];
    vs[3] = parent.vertices[3];

    vs[4] = (vs[0] + vs[2]) / 2.0;
    vs[5] = (vs[1] + vs[3]) / 2.0;
    vs[6] = (vs[0] + vs[1]) / 2.0;
    vs[7] = (vs[2] + vs[3]) / 2.0;

    vs[8] = (vs[0] + vs[1] + vs[2] + vs[3]) / 4.0;

    Quad vl0 = {vs[0], vs[6], vs[4], vs[8]};
    Quad vl1 = {vs[6], vs[1], vs[8], vs[5]};
    Quad vl2 = {vs[4], vs[8], vs[2], vs[7]};
    Quad vl3 = {vs[8], vs[5], vs[7], vs[3]};

    return {AmrCell(vl0), AmrCell(vl1), AmrCell(vl2), AmrCell(vl3)};
}

/// @brief Создать геометрию дочерних ячеек по родительским вершинам (2D)
/// @param parent_vertices Вершины родительской ячейки
/// @return Массив с дочерними ячейками
template <>
std::array<geom::AmrCell, CpC(2)> create_children<2>(AmrCell& parent) {
    using geom::SqQuad;
    using geom::AmrCell;

    const AmrVertices &vertices = parent.vertices;

    // Собираем отображение ячейки
    SqQuad quad = {
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

    // @formatter:off
    SqQuad vl0 = {
            quad.vs<-1, -1>(), quad(-0.5, -1.0), quad.vs<0, -1>(),
            quad(-1.0, -0.5),  quad(-0.5, -0.5), quad(0.0, -0.5),
            quad.vs<-1,  0>(), quad(-0.5,  0.0), quad.vs<0,  0>()
    };
    SqQuad vl1 = {
            quad.vs<0, -1>(), quad(0.5, -1.0), quad.vs<1, -1>(),
            quad(0.0, -0.5),  quad(0.5, -0.5), quad(1.0, -0.5),
            quad.vs<0,  0>(), quad(0.5,  0.0), quad.vs<1,  0>()
    };
    SqQuad vl2 = {
            quad.vs<-1, 0>(), quad(-0.5, 0.0), quad.vs<0, 0>(),
            quad(-1.0, 0.5),  quad(-0.5, 0.5), quad(0.0, 0.5),
            quad.vs<-1, 1>(), quad(-0.5, 1.0), quad.vs<0, 1>()
    };
    SqQuad vl3 = {
            quad.vs<0, 0>(), quad(0.5, 0.0), quad.vs<1, 0>(),
            quad(0.0, 0.5),  quad(0.5, 0.5), quad(1.0, 0.5),
            quad.vs<0, 1>(), quad(0.5, 1.0), quad.vs<1, 1>()
    };
    // @formatter:on

    return {AmrCell(vl0), AmrCell(vl1), AmrCell(vl2), AmrCell(vl3)};
}

/// @brief Создать геометрию дочернх ячеек по родительским вершинам (3D)
/// @param parent_vertices Вершины родительской ячейки
/// @return Массив с дочерними ячейками
template <>
std::array<AmrCell, CpC(3)> create_children<3>(AmrCell& parent) {
    SqCube& vs = (SqCube&) parent.vertices;
   
    Cube vl1 = {vs[iww(0, 0, 0)], vs[iww(1, 0, 0)], vs[iww(0, 1, 0)], vs[iww(1, 1, 0)],
                vs[iww(0, 0, 1)], vs[iww(1, 0, 1)], vs[iww(0, 1, 1)], vs[iww(1, 1, 1)]};
    Cube vl2 = {vs[iww(1, 0, 0)], vs[iww(2, 0, 0)], vs[iww(1, 1, 0)], vs[iww(2, 1, 0)],
                vs[iww(1, 0, 1)], vs[iww(2, 0, 1)], vs[iww(1, 1, 1)], vs[iww(2, 1, 1)]};
    Cube vl3 = {vs[iww(0, 1, 0)], vs[iww(1, 1, 0)], vs[iww(0, 2, 0)], vs[iww(1, 2, 0)],
                vs[iww(0, 1, 1)], vs[iww(1, 1, 1)], vs[iww(0, 2, 1)], vs[iww(1, 2, 1)]};
    Cube vl4 = {vs[iww(1, 1, 0)], vs[iww(2, 1, 0)], vs[iww(1, 2, 0)], vs[iww(2, 2, 0)],
                vs[iww(1, 1, 1)], vs[iww(2, 1, 1)], vs[iww(1, 2, 1)], vs[iww(2, 2, 1)]};
    Cube vl5 = {vs[iww(0, 0, 1)], vs[iww(1, 0, 1)], vs[iww(0, 1, 1)], vs[iww(1, 1, 1)],
                vs[iww(0, 0, 2)], vs[iww(1, 0, 2)], vs[iww(0, 1, 2)], vs[iww(1, 1, 2)]};
    Cube vl6 = {vs[iww(1, 0, 1)], vs[iww(2, 0, 1)], vs[iww(1, 1, 1)], vs[iww(2, 1, 1)],
                vs[iww(1, 0, 2)], vs[iww(2, 0, 2)], vs[iww(1, 1, 2)], vs[iww(2, 1, 2)]};
    Cube vl7 = {vs[iww(0, 1, 1)], vs[iww(1, 1, 1)], vs[iww(0, 2, 1)], vs[iww(1, 2, 1)],
                vs[iww(0, 1, 2)], vs[iww(1, 1, 2)], vs[iww(0, 2, 2)], vs[iww(1, 2, 2)]};
    Cube vl8 = {vs[iww(1, 1, 1)], vs[iww(2, 1, 1)], vs[iww(1, 2, 1)], vs[iww(2, 2, 1)],
                vs[iww(1, 1, 2)], vs[iww(2, 1, 2)], vs[iww(1, 2, 2)], vs[iww(2, 2, 2)]};

    Vector3d zero(0.0, 0.0, 0.0);
    Cube vl0 = {zero, zero, zero, zero, zero, zero, zero, zero};

    return {
            AmrCell(vl1), AmrCell(vl2), AmrCell(vl3), AmrCell(vl4),
            AmrCell(vl5), AmrCell(vl6), AmrCell(vl7), AmrCell(vl8)
    };
}

/// @brief Создать дочерние ячейки
/// @param cell Родительская ячейка
/// @param ic Индекс родительской ячейки
/// @return Массив с дочерними ячейками, дочерние ячейки имеют законченный вид
/// (необходимое число граней, правильную линковку (на старные ячейки))
template<int dim>
std::array<AmrCell, CpC(dim)> get_children(AmrCell &cell) {
    auto children = create_children<dim>(cell);

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
                AmrFace &child_face = child.faces[side];
                auto child_fc = child_face.center<dim>(child.vertices);

                for (auto s: subface_sides<dim>(side)) {
                    AmrFace& cell_face = cell.faces[s];
                    auto cell_fc = cell_face.center<dim>(cell.vertices);

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
std::array<Storage::Item, CpC(dim)> select_children(
        Storage& cells, const std::array<AmrCell, CpC(dim)>& children);

/// @brief Возвращает итераторы дочерних ячеек (2D)
template <>
std::array<Storage::Item, CpC(2)> select_children<2>(
        Storage& cells, const std::array<AmrCell, CpC(2)>& children) {
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
        Storage& cells, const std::array<AmrCell, CpC(3)>& children) {
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
    auto item = locals[ic];
    AmrCell& parent = item.geom();

    auto children = get_children<dim>(parent);

    for (int i = 0; i < CpC(dim); ++i) {
        AmrCell& child = locals[parent.next + i];

        // Возможна оптимизация (!), создавать дочерние ячейки сразу
        // на нужном месте в хранилище
        children[i].copy_to(child);

        for (int s = 0; s < FpC(dim); ++s) {
            AmrFace &f1 = child.faces[s];
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

            if (nei_wanted_lvl > child.level) {
                split_face<dim>(child, Side(s));
            }
        }
    }

    op.split<dim>(item, select_children<dim>(locals, children));

    parent.set_undefined();
}

} // namespace amr
} // namespace mesh
} // namespace zephyr