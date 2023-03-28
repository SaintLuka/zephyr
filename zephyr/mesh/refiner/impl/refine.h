/// @file Файл содержит реализацию функций для разбиения ячейки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/math/geom/cell.h>
#include <zephyr/math/geom/maps.h>

#include <zephyr/mesh/refiner/impl/common.h>
#include <zephyr/mesh/refiner/impl/faces.h>
#include <zephyr/mesh/refiner/impl/coarse.h>

namespace zephyr { namespace mesh { namespace impl {

using namespace zephyr::data;
using namespace zephyr::math;
using amrData;

/// @brief Создать геометрию дочернх ячеек по родительской ячейке
/// @param parent Родительская ячейка
/// @return Массив с дочерними ячейками
template <unsigned int dim>
std::array<geom::Cell, CpC(dim)> create_children(Storage::iterator parent);

/// @brief Создать геометрию дочерних ячеек по родительской ячейке
/// @param parent Родительская ячейка
/// @return Массив с дочерними ячейками
std::array<geom::Cell, CpC(2)> create_children_simple(Storage::iterator parent) {
    using geom::ShortList2D;
    using geom::Cell;

    std::array<Vector3d, 9> vs;
    vs[0] = (Vector3d &) parent[vertices].list[0];
    vs[1] = (Vector3d &) parent[vertices].list[1];
    vs[2] = (Vector3d &) parent[vertices].list[2];
    vs[3] = (Vector3d &) parent[vertices].list[3];

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
std::array<geom::Cell, CpC(2)> create_children<2>(Storage::iterator parent) {
    using geom::LargeList2D;
    using geom::Mapping2D;
    using geom::Cell;

    const _vertices_& vertices = parent[vertices];

    // Собираем отображение ячейки
    LargeList2D vs = {
            (Vector3d &) vertices.list[0],
            (Vector3d &) vertices.list[1],
            (Vector3d &) vertices.list[2],
            (Vector3d &) vertices.list[3],
            (Vector3d &) vertices.list[4],
            (Vector3d &) vertices.list[5],
            (Vector3d &) vertices.list[6],
            (Vector3d &) vertices.list[7],
            (Vector3d &) vertices.list[8]
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
std::array<geom::Cell, CpC(3)> create_children<3>(Storage::iterator parent) {
    using geom::ShortList3D;
    using geom::LargeList3D;
    using topology::iww;
    using geom::Cell;

    auto& vertices = parent[vertices];
    
    LargeList3D& vs = (LargeList3D&) parent[vertices].list;
   
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
template<unsigned int dim>
std::array<geom::Cell, CpC(dim)> get_children(Storage::iterator &cell, size_t ic, unsigned int rank) {
    using side;
    using Vector3d;

    auto children = create_children<dim>(cell);

    auto children_by_side = get_children_by_side<dim>();

    // По умолчанию дети ссылаются на родительскую ячейку
    for (auto &child: children) {
        for (unsigned int s = 0; s < FpC(dim); ++s) {
            child.faces[s].adjacent.rank = rank;
            child.faces[s].adjacent.index = ic;
            child.faces[s].adjacent.ghost = std::numeric_limits<unsigned int>::max();
        }
    }

    for (unsigned int i = 0; i < CpC(dim); ++i) {
        children[i].amrData.next = cell[amrData].next + i;
        children[i].amrData.flag = 0;
        children[i].amrData.base_id = cell[amrData].base_id;
        children[i].amrData.level = cell[amrData].level + 1;
        children[i].amrData.z = CpC(dim) * cell[amrData].z + i;
    }

    // Выставить граничный флаг
    for (unsigned int side = 0; side < FpC(dim); ++side) {
        auto flag = cell[faces].list[side].boundary;
        for (int i: children_by_side[side]) {
            children[i].faces[side].boundary = flag;
        }
    }

    // Далее необходимо слинковать дочерние ячейки с соседями
    for (unsigned int side = 0; side < FpC(dim); ++side) {
        auto flag = cell[faces].list[side].boundary;
        for (int i: children_by_side[side]) {
            children[i].faces[side].boundary = flag;
        }

        if (cell[faces].list[side + 6].is_undefined()) {
            // Ячейка имела простую грань
            auto adj = cell[faces].list[side].adjacent;
            for (int i: children_by_side[side]) {
                children[i].faces[side].adjacent = adj;
            }
        } else {
            // Ячейка имела сложную грань
            for (int i: children_by_side[side]) {
                geom::Cell &child = children[i];
                _face_ &child_face = child.faces[side];
                auto child_fc = face_center<dim>(child_face, child.vertices);

                for (auto s: subface_sides<dim>(side)) {
                    _face_ &cell_face = cell[faces].list[s];
                    auto cell_fc = face_center<dim>(cell_face, cell[vertices]);

                    if (distance(child_fc, cell_fc) < 1.0e-5 * cell[size]) {
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
template <unsigned int dim>
std::array<Storage::iterator, CpC(dim)> select_children(
        Storage& cells, const std::array<geom::Cell, CpC(dim)>& children);

/// @brief Возвращает итераторы дочерних ячеек (2D)
template <>
std::array<Storage::iterator, CpC(2)> select_children<2>(
        Storage& cells, const std::array<geom::Cell, CpC(2)>& children) {
    return {
            cells[children[0].amrData.next],
            cells[children[1].amrData.next],
            cells[children[2].amrData.next],
            cells[children[3].amrData.next]
    };
}

/// @brief Возвращает итераторы дочерних ячеек (3D)
template <>
std::array<Storage::iterator, CpC(3)> select_children<3>(
        Storage& cells, const std::array<geom::Cell, CpC(3)>& children) {
    return {
            cells[children[0].amrData.next],
            cells[children[1].amrData.next],
            cells[children[2].amrData.next],
            cells[children[3].amrData.next],
            cells[children[4].amrData.next],
            cells[children[5].amrData.next],
            cells[children[6].amrData.next],
            cells[children[7].amrData.next]
    };
}

/// @brief Производит разбиение ячейки, дочерние ячейки помещает в Storage
/// @param locals Хранилище ячеек
/// @param ic Индекс родительской ячейки в хранилище
/// @param op Оператор разделения данных
/// @details Дочерние ячейки правильно ссылаются друг на друга, на гранях
/// adjacent указан на старые ячейки.
template<unsigned int dim>
void refine_cell(Storage &locals, Storage &aliens, unsigned int rank, size_t ic, const DataDistributor& op) {
    auto cell = locals[ic];

    cell[element].set_undefined();

    auto children = get_children<dim>(cell, ic, rank);

    for (unsigned int i = 0; i < CpC(dim); ++i) {
        auto j = locals[ic][amrData].next + i;

        copy_data(cell, locals[j]);

        locals[j][coords] = children[i].coords;
        locals[j][vertices] = children[i].vertices;
        locals[j][faces] = children[i].faces;
        locals[j][size] = children[i].size;
        locals[j][amrData] = children[i].amrData;
        locals[j][element].kind = kind::EULER;
        locals[j][element].dimension = dim;
        locals[j][element].index = std::numeric_limits<size_t>::max();
        locals[j][element].rank = rank;

        for (unsigned int s = 0; s < FpC(dim); ++s) {
            _face_ &f1 = locals[j][faces].list[s];
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
                nei_wanted_lvl = neib[amrData].level + neib[amrData].flag;
            }
            else {
                // Удаленная ячейка
#if SCRUTINY
                if (adj.ghost >= aliens.size()) {
                    throw std::runtime_error("adjacent.ghost out of range (refine_cell)");
                }
#endif
                auto neib = aliens[adj.ghost];
                nei_wanted_lvl = neib[amrData].level + neib[amrData].flag;
            }

            if (nei_wanted_lvl > locals[j][amrData].level) {
                split_face<dim>(locals[j], side(s));
            }
        }
    }

    op.split<dim>(cell, select_children<dim>(locals, children));
}

} // namespace impl
} // namespace mesh
} // namespace zephyr