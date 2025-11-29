// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#include <zephyr/mesh/amr/common.h>

namespace zephyr::mesh::amr {

/// @brief Устанавливает на новых подгранях индексы вершин, копирует значения
/// boundary и adjacent с исходной грани.
/// @param faces Массив граней
/// @param face_beg Индекс первой грани ячейки
template <int dim, int int_side>
void setup_faces_topo(AmrFaces& faces, index_t face_beg) {
    constexpr Side<dim> side{int_side};

    // Индекс исходной грани
    const index_t orig_face = face_beg + side;

    // Исходная грань превращается в первую подгрань
    faces.vertices[orig_face] = side[0].cf();

    const index_t subface_1 = face_beg + side[1];
    faces.vertices[subface_1] = side[1].cf();
    faces.boundary[subface_1] = faces.boundary[orig_face];
    faces.adjacent.rank [subface_1] = faces.adjacent.rank [orig_face];
    faces.adjacent.index[subface_1] = faces.adjacent.index[orig_face];
    faces.adjacent.alien[subface_1] = faces.adjacent.alien[orig_face];
    faces.adjacent.basic[subface_1] = faces.adjacent.basic[orig_face];

    if constexpr (dim > 2) {
        const index_t subface_2 = face_beg + side[2];
        faces.vertices[subface_2] = side[2].cf();
        faces.boundary[subface_2] = faces.boundary[orig_face];
        faces.adjacent.rank [subface_2] = faces.adjacent.rank [orig_face];
        faces.adjacent.index[subface_2] = faces.adjacent.index[orig_face];
        faces.adjacent.alien[subface_2] = faces.adjacent.alien[orig_face];
        faces.adjacent.basic[subface_2] = faces.adjacent.basic[orig_face];

        const index_t subface_3 = face_beg + side[3];
        faces.vertices[subface_3] = side[3].cf();
        faces.boundary[subface_3] = faces.boundary[orig_face];
        faces.adjacent.rank [subface_3] = faces.adjacent.rank [orig_face];
        faces.adjacent.index[subface_3] = faces.adjacent.index[orig_face];
        faces.adjacent.alien[subface_3] = faces.adjacent.alien[orig_face];
        faces.adjacent.basic[subface_3] = faces.adjacent.basic[orig_face];
    }
}

/// @brief Устанавливает нормали, площадь и центр новой подграни
/// @param faces Массив граней
/// @param iface Индекс грани, для которой необходимо посчитать геметрию
/// @param vertices Вершины ячейки
template <int dim, bool axial>
void setup_face_geom(AmrFaces& faces, index_t iface, const SqMap<dim>& vertices) {
    // Точка внутри ячейки
    const auto& C = vertices.center();

    if constexpr (dim == 2) {
        // Вершины грани
        Line vl = {
            vertices[faces.vertices[iface][0]],
            vertices[faces.vertices[iface][1]],
        };

        // Установить длину и нормаль
        faces.area[iface]   = vl.length();
        if constexpr (!axial) {
            faces.center[iface] = vl.centroid();
        }
        else {
            faces.center[iface] = vl.centroid(axial);
        }
        faces.normal[iface] = vl.normal(C);

        if constexpr (axial) {
            faces.area_alt[iface] = vl.area_as();
        }
    }
    else {
        // Вершины подграней
        Quad vl = {
            vertices[faces.vertices[iface][0]],
            vertices[faces.vertices[iface][1]],
            vertices[faces.vertices[iface][2]],
            vertices[faces.vertices[iface][3]],
        };

        // Установить площадь и нормаль
        faces.area[iface]   = vl.area();
        faces.center[iface] = vl.center();
        faces.normal[iface] = vl.normal(C);
    }
}

/// @brief Устанавливает геометрию для всех подграней
/// @param faces Массив граней
/// @param face_beg Индекс первой грани ячейки
/// @param vertices Вершины ячейки
template <int dim, int side_in, bool axial>
void setup_faces_geom(AmrFaces& faces, index_t face_beg, const SqMap<dim>& vertices) {
    constexpr Side<dim> side{side_in};

    setup_face_geom<dim, axial>(faces, face_beg + side[0], vertices);
    setup_face_geom<dim, axial>(faces, face_beg + side[1], vertices);

    if constexpr (dim > 2) {
        setup_face_geom<dim, false>(faces, face_beg + side[2], vertices);
        setup_face_geom<dim, false>(faces, face_beg + side[3], vertices);
    }
}

/// @brief Разбить простую грань на две или четыре подграни. Устанавливает
/// геометрию подграней, индексы смежности копирует с исходной грани.
/// @param faces Массив граней
/// @param face_beg Индекс первой грани ячейки
/// @param vertices Вершины ячейки
template <int dim, int side, bool axial = false>
void split_face_impl(AmrFaces& faces, index_t face_beg, const SqMap<dim>& vertices) {
#if SCRUTINY
    if (faces.is_actual(face_beg + Side<dim>(side)[1])) {
        throw std::runtime_error("Attempt to split complex face");
    }
#endif
    setup_faces_topo<dim, side>(faces, face_beg);
    setup_faces_geom<dim, side, axial>(faces, face_beg, vertices);
}

/// @brief Объединяет подграни в одну грань. Устанавливает геометрию новой грани,
/// подграни делает неактуальными, индексы смежности остаются с первой подграни.
/// @param faces Массив граней
/// @param face_beg Индекс первой грани ячейки
/// @param vertices Вершины ячейки
template<int dim, int int_side, bool axial>
void merge_faces_impl(AmrFaces& faces, index_t face_beg, const SqMap<dim>& vertices) {
    constexpr Side<dim> side{int_side};
#if SCRUTINY
    if (faces.is_undefined(face_beg + side[1])) {
        throw std::runtime_error("Attempt to split complex face");
    }
#endif

    faces.vertices[face_beg + side] = side.sf();

    faces.set_undefined(face_beg + side[1]);
    if constexpr (dim > 2) {
        faces.set_undefined(face_beg + side[2]);
        faces.set_undefined(face_beg + side[3]);
    }
    setup_face_geom<dim, axial>(faces, face_beg + side, vertices);
}

/// @brief Разбить простую грань на две или четыре подграни. Устанавливает
/// геометрию подграней, индексы смежности копирует с исходной грани.
/// @param cells Хранилище ячеек
/// @param ic Индекс целевой ячейки в хранилище
/// @param side Простая грань ячейки
template <int dim>
void split_face(AmrCells& cells, index_t ic, Side<dim> side) {
    scrutiny_check(side < Side<dim>::count(), "Attempt to split subface");

    using split_face_2D_t = void (*)(AmrFaces&, index_t, const SqQuad&);
    static constexpr split_face_2D_t split_lookup_table_2D[4][2] = {
        {split_face_impl<2, Side2D::L, false>, split_face_impl<2, Side2D::L, true>},
        {split_face_impl<2, Side2D::R, false>, split_face_impl<2, Side2D::R, true>},
        {split_face_impl<2, Side2D::B, false>, split_face_impl<2, Side2D::B, true>},
        {split_face_impl<2, Side2D::T, false>, split_face_impl<2, Side2D::T, true>},
    };

    using split_face_3D_t = void (*)(AmrFaces&, index_t, const SqCube&);
    static constexpr split_face_3D_t split_lookup_table_3D[6] = {
        split_face_impl<3, Side3D::L, false>,
        split_face_impl<3, Side3D::R, false>,
        split_face_impl<3, Side3D::B, false>,
        split_face_impl<3, Side3D::T, false>,
        split_face_impl<3, Side3D::Z, false>,
        split_face_impl<3, Side3D::F, false>
    };

    if constexpr (dim == 2) {
        return split_lookup_table_2D[side][cells.axial()](
            cells.faces, cells.face_begin[ic], cells.mapping<dim>(ic));
    }
    else {
        return split_lookup_table_3D[side](
            cells.faces, cells.face_begin[ic], cells.mapping<dim>(ic));
    }
}

/// @brief Объединяет подграни в одну грань. Устанавливает геометрию новой грани,
/// подграни делает неактуальными, индексы смежности остаются с первой подграни.
/// @param cells Хранилище ячеек
/// @param ic Индекс целевой ячейки в хранилище
/// @param side Простая грань ячейки
template<int dim>
void merge_faces(AmrCells& cells, index_t ic, Side<dim> side) {
    scrutiny_check(side < Side<dim>::count(), "Attempt to merge subface");

    using merge_faces_2D_t = void (*)(AmrFaces&, index_t, const SqQuad&);
    static constexpr merge_faces_2D_t merge_lookup_table_2D[4][2] = {
        {merge_faces_impl<2, Side2D::L, false>, merge_faces_impl<2, Side2D::L, true>},
        {merge_faces_impl<2, Side2D::R, false>, merge_faces_impl<2, Side2D::R, true>},
        {merge_faces_impl<2, Side2D::B, false>, merge_faces_impl<2, Side2D::B, true>},
        {merge_faces_impl<2, Side2D::T, false>, merge_faces_impl<2, Side2D::T, true>},
    };

    using merge_faces_3D_t = void (*)(AmrFaces&, index_t, const SqCube&);
    static constexpr merge_faces_3D_t merge_lookup_table_3D[6] = {
        merge_faces_impl<3, Side3D::L, false>,
        merge_faces_impl<3, Side3D::R, false>,
        merge_faces_impl<3, Side3D::B, false>,
        merge_faces_impl<3, Side3D::T, false>,
        merge_faces_impl<3, Side3D::Z, false>,
        merge_faces_impl<3, Side3D::F, false>
    };

    if constexpr (dim == 2) {
        return merge_lookup_table_2D[side][cells.axial()](
            cells.faces, cells.face_begin[ic], cells.mapping<dim>(ic));
    }
    else {
        return merge_lookup_table_3D[side](
            cells.faces, cells.face_begin[ic], cells.mapping<dim>(ic));
    }    
}

} // namespace zephyr::mesh::amr