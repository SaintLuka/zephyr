/// @file Файл содержит операции разбиения и огрубления для граней ячейки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/mesh/primitives/Side3D.h>
#include <zephyr/mesh/euler/soa_mesh.h>
#include <zephyr/mesh/amr2/common.h>


namespace zephyr::mesh::amr2 {

template <int dim>
using SqMap = std::conditional_t<dim < 3, SqQuad, SqCube>;

/// @brief Устанавливает на новых подгранях значение adjacent и
/// boundary со старой грани
template <int dim, int side>
void split_face_prepare(index_t iface, SoaCell::SoaFace& faces) {
    faces.boundary[iface + (side + 6)] = faces.boundary[side];
    faces.adjacent.rank[iface + (side + 6)] = faces.adjacent.rank[side];
    faces.adjacent.local_index[iface + (side + 6)] = faces.adjacent.local_index[side];
    faces.adjacent.owner_index[iface + (side + 6)] = faces.adjacent.owner_index[side];

    if constexpr (dim > 2) {
        faces.boundary[iface + (side + 12)] = faces.boundary[side];
        faces.adjacent.rank[iface + (side + 12)] = faces.adjacent.rank[side];
        faces.adjacent.local_index[iface + (side + 12)] = faces.adjacent.local_index[side];
        faces.adjacent.owner_index[iface + (side + 12)] = faces.adjacent.owner_index[side];

        faces.boundary[iface + (side + 18)] = faces.boundary[side];
        faces.adjacent.rank[iface + (side + 18)] = faces.adjacent.rank[side];
        faces.adjacent.local_index[iface + (side + 18)] = faces.adjacent.local_index[side];
        faces.adjacent.owner_index[iface + (side + 18)] = faces.adjacent.owner_index[side];
    }
}

/// @brief Устанавливает на новых подгранях индексы вершин.
template <int dim, int side>
void split_face_indices(index_t iface, SoaCell::SoaFace& faces) {
    faces.vertices[iface + side] = face_indices::cf<dim, Side3D(side)>();
    faces.vertices[iface + (side + 6)] = face_indices::cf<dim, Side3D(side + 6)>();

    if constexpr (dim > 2) {
        faces.vertices[iface + (side + 12)] = face_indices::cf<dim, Side3D(side + 12)>();
        faces.vertices[iface + (side + 18)] = face_indices::cf<dim, Side3D(side + 18)>();
    }
}

/// @brief Устанавливает нормали и площади новых подграней
template <int dim>
void setup_face_features(index_t iface, SoaCell::SoaFace& faces, const SqMap<dim>& vertices, bool axial = false);

template <>
inline void setup_face_features<2>(index_t iface, SoaCell::SoaFace& faces, const SqQuad& vertices, bool axial) {
    // Точка внутри ячейки
    Vector3d C = vertices.vs<0, 0>();

    // Вершины грани
    Line vl = {
            vertices[faces.vertices[iface][0]],
            vertices[faces.vertices[iface][1]],
    };

    // Установить длину и нормаль
    faces.area[iface]   = vl.length();
    faces.center[iface] = vl.centroid(axial);
    faces.normal[iface] = vl.normal(C);

    if (axial) { faces.area_alt[iface] = vl.area_as(); }
}

template <>
inline void setup_face_features<3>(index_t iface, SoaCell::SoaFace& faces, const SqCube& vertices, bool axial) {
    // Точка внутри ячейки
    const Vector3d& C = vertices.vs<0, 0, 0>();

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


/// @brief Устанавливает нормали и площади новых подграней
template <int dim>
void split_face_features(index_t iface, SoaCell::SoaFace& faces,
    const SqMap<dim>& vertices, int side, bool axial) {
    setup_face_features<dim>(iface + side, faces, vertices, axial);
    setup_face_features<dim>(iface + side + 6, faces, vertices, axial);

    if constexpr (dim > 2) {
        setup_face_features<dim>(iface + side, faces, vertices);
        setup_face_features<dim>(iface + side, faces, vertices);
    }
}

/// @brief Разбивает грань на стороне side на подграни.
/// Устанавливает площади и нормали новых граней.
/// Размерность и сторона являются аргументами шаблона
template <int dim, int side>
void split_face(index_t iface, SoaCell::SoaFace& faces, const SqMap<dim>& vertices, bool axial = false) {
#if SCRUTINY
    if (faces.is_actual(iface + side + 6)) {
        throw std::runtime_error("Try to cut complex face");
    }
#endif

    split_face_prepare<dim, side>(iface, faces);

    split_face_indices<dim, side>(iface, faces);

    split_face_features<dim>(iface, faces, vertices, side, axial);
}

/// @brief Разбивает грань на стороне side на подграни.
/// Устанавливает площади и нормали новых граней.
/// Преобразует аргумент функции side в аргумент шаблона.
template <int dim>
void split_face(index_t iface, SoaCell::SoaFace& faces, const SqMap<dim>& vertices, int side, bool axial);

template <>
inline void split_face<2>(index_t iface, SoaCell::SoaFace& faces, const SqQuad& vertices, int side, bool axial) {
    switch (side) {
        case Side3D::LEFT:
            split_face<2, Side3D::L>(iface, faces, vertices, axial);
            break;
        case Side3D::RIGHT:
            split_face<2, Side3D::R>(iface, faces, vertices, axial);
            break;
        case Side3D::BOTTOM:
            split_face<2, Side3D::B>(iface, faces, vertices, axial);
            break;
        default:
            split_face<2, Side3D::T>(iface, faces, vertices, axial);
            break;
    }
}

template <>
inline void split_face<3>(index_t iface, SoaCell::SoaFace& faces, const SqCube& vertices, int side, bool axial) {
    switch (side) {
        case Side3D::LEFT:
            split_face<3, Side3D::L>(iface, faces, vertices);
            break;
        case Side3D::RIGHT:
            split_face<3, Side3D::R>(iface, faces, vertices);
            break;
        case Side3D::BOTTOM:
            split_face<3, Side3D::B>(iface, faces, vertices);
            break;
        case Side3D::TOP:
            split_face<3, Side3D::T>(iface, faces, vertices);
            break;
        case Side3D::BACK:
            split_face<3, Side3D::X>(iface, faces, vertices);
            break;
        default:
            split_face<3, Side3D::F>(iface, faces, vertices);
            break;
    }
}

/// @brief Разбивает грань на стороне side на подграни.
/// Устанавливает площади и нормали новых граней.
template <int dim>
void split_face(index_t ic, SoaCell& cells, int side) {
    split_face<dim>(cells.face_begin[ic] + side, cells.faces, cells.get_vertices<dim>(ic), side, cells.axial);
}

/// @brief Объединяет подграни в одну грань на стороне side
template<int dim, int side>
void merge_faces(index_t ic, SoaCell& cells) {
    index_t iface  = cells.face_begin[ic] + side;
    const SqMap<dim> &vertices = cells.get_vertices<dim>(ic);

    cells.faces.vertices[iface + side] = face_indices::sf<dim, side>();

    cells.faces.set_undefined(iface + (side + 6));

    if constexpr (dim > 2) {
        cells.faces.set_undefined(iface + (side + 12));
        cells.faces.set_undefined(iface + (side + 18));
    }

    setup_face_features<dim>(iface + side, cells.faces, vertices, cells.axial);
}

/// @brief Объединяет подграни в одну грань на стороне side
template<int dim>
void merge_faces(index_t ic, SoaCell& cells, Side3D side);

template<>
inline void merge_faces<2>(index_t ic, SoaCell& cells, Side3D side) {
    switch (side) {
        case Side3D::LEFT:
            merge_faces<2, Side3D::L>(ic, cells);
            break;
        case Side3D::RIGHT:
            merge_faces<2, Side3D::R>(ic, cells);
            break;
        case Side3D::BOTTOM:
            merge_faces<2, Side3D::B>(ic, cells);
            break;
        default:
            merge_faces<2, Side3D::T>(ic, cells);
            break;
    }
}

template<>
inline void merge_faces<3>(index_t ic, SoaCell& cells, Side3D side) {
    switch (side) {
        case Side3D::LEFT:
            merge_faces<3, Side3D::L>(ic, cells);
            break;
        case Side3D::RIGHT:
            merge_faces<3, Side3D::R>(ic, cells);
            break;
        case Side3D::BOTTOM:
            merge_faces<3, Side3D::B>(ic, cells);
            break;
        case Side3D::TOP:
            merge_faces<3, Side3D::T>(ic, cells);
            break;
        case Side3D::BACK:
            merge_faces<3, Side3D::X>(ic, cells);
            break;
        default:
            merge_faces<3, Side3D::F>(ic, cells);
            break;
    }
}

} // namespace zephyr::mesh::amr2