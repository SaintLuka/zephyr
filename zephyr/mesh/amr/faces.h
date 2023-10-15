/// @file Файл содержит операции разбиения и огрубления для граней ячейки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/mesh/amr/common.h>

namespace zephyr { namespace mesh { namespace amr {

/// @brief Устанавливает на новых подгранях значение adjacent и
/// boundary со старой грани
template <int dim, int side>
void split_face_prepare(AmrFaces& faces) {
    faces[side + 6].boundary = faces[side].boundary;
    faces[side + 6].adjacent = faces[side].adjacent;

    if (dim > 2) {
        faces[side + 12].boundary = faces[side].boundary;
        faces[side + 12].adjacent = faces[side].adjacent;
        faces[side + 18].boundary = faces[side].boundary;
        faces[side + 18].adjacent = faces[side].adjacent;
    }
}

/// @brief Устанавливает на новых подгранях индексы вершин.
template <int dim, int side>
void split_face_indices(AmrFaces& faces);

template <>
void split_face_indices<2, Side::LEFT>(AmrFaces& faces) {
    faces[Side::LEFT0].vertices[0] = iww(0, 0);
    faces[Side::LEFT0].vertices[1] = iww(0, 1);
    faces[Side::LEFT1].vertices[0] = iww(0, 1);
    faces[Side::LEFT1].vertices[1] = iww(0, 2);
}

template <>
void split_face_indices<2, Side::RIGHT>(AmrFaces& faces) {
    faces[Side::RIGHT0].vertices[0] = iww(2, 0);
    faces[Side::RIGHT0].vertices[1] = iww(2, 1);
    faces[Side::RIGHT1].vertices[0] = iww(2, 1);
    faces[Side::RIGHT1].vertices[1] = iww(2, 2);
}

template <>
void split_face_indices<2, Side::BOTTOM>(AmrFaces& faces) {
    faces[Side::BOTTOM0].vertices[0] = iww(0, 0);
    faces[Side::BOTTOM0].vertices[1] = iww(1, 0);
    faces[Side::BOTTOM1].vertices[0] = iww(1, 0);
    faces[Side::BOTTOM1].vertices[1] = iww(2, 0);
}

template <>
void split_face_indices<2, Side::TOP>(AmrFaces& faces) {
    faces[Side::TOP0].vertices[0] = iww(0, 2);
    faces[Side::TOP0].vertices[1] = iww(1, 2);
    faces[Side::TOP1].vertices[0] = iww(1, 2);
    faces[Side::TOP1].vertices[1] = iww(2, 2);
}

template <>
void split_face_indices<3, Side::LEFT>(AmrFaces& faces) {
    faces[Side::LEFT0].vertices = {iww(0, 0, 0), iww(0, 1, 0), iww(0, 0, 1), iww(0, 1, 1)};
    faces[Side::LEFT1].vertices = {iww(0, 1, 0), iww(0, 2, 0), iww(0, 1, 1), iww(0, 2, 1)};
    faces[Side::LEFT2].vertices = {iww(0, 0, 1), iww(0, 1, 1), iww(0, 0, 2), iww(0, 1, 2)};
    faces[Side::LEFT3].vertices = {iww(0, 1, 1), iww(0, 2, 1), iww(0, 1, 2), iww(0, 2, 2)};
}

template <>
void split_face_indices<3, Side::RIGHT>(AmrFaces& faces) {
    faces[Side::RIGHT0].vertices = {iww(2, 0, 0), iww(2, 1, 0), iww(2, 0, 1), iww(2, 1, 1)};
    faces[Side::RIGHT1].vertices = {iww(2, 1, 0), iww(2, 2, 0), iww(2, 1, 1), iww(2, 2, 1)};
    faces[Side::RIGHT2].vertices = {iww(2, 0, 1), iww(2, 1, 1), iww(2, 0, 2), iww(2, 1, 2)};
    faces[Side::RIGHT3].vertices = {iww(2, 1, 1), iww(2, 2, 1), iww(2, 1, 2), iww(2, 2, 2)};
}

template <>
void split_face_indices<3, Side::BOTTOM>(AmrFaces& faces) {
    faces[Side::BOTTOM0].vertices = {iww(0, 0, 0), iww(1, 0, 0), iww(0, 0, 1), iww(1, 0, 1)};
    faces[Side::BOTTOM1].vertices = {iww(1, 0, 0), iww(2, 0, 0), iww(1, 0, 1), iww(2, 0, 1)};
    faces[Side::BOTTOM2].vertices = {iww(0, 0, 1), iww(1, 0, 1), iww(0, 0, 2), iww(1, 0, 2)};
    faces[Side::BOTTOM3].vertices = {iww(1, 0, 1), iww(2, 0, 1), iww(1, 0, 2), iww(2, 0, 2)};
}

template <>
void split_face_indices<3, Side::TOP>(AmrFaces& faces) {
    faces[Side::TOP0].vertices = {iww(0, 2, 0), iww(1, 2, 0), iww(0, 2, 1), iww(1, 2, 1)};
    faces[Side::TOP1].vertices = {iww(1, 2, 0), iww(2, 2, 0), iww(1, 2, 1), iww(2, 2, 1)};
    faces[Side::TOP2].vertices = {iww(0, 2, 1), iww(1, 2, 1), iww(0, 2, 2), iww(1, 2, 2)};
    faces[Side::TOP3].vertices = {iww(1, 2, 1), iww(2, 2, 1), iww(1, 2, 2), iww(2, 2, 2)};
}

template <>
void split_face_indices<3, Side::BACK>(AmrFaces& faces) {
    faces[Side::BACK0].vertices = {iww(0, 0, 0), iww(1, 0, 0), iww(0, 1, 0), iww(1, 1, 0)};
    faces[Side::BACK1].vertices = {iww(1, 0, 0), iww(2, 0, 0), iww(1, 1, 0), iww(2, 1, 0)};
    faces[Side::BACK2].vertices = {iww(0, 1, 0), iww(1, 1, 0), iww(0, 2, 0), iww(1, 2, 0)};
    faces[Side::BACK3].vertices = {iww(1, 1, 0), iww(2, 1, 0), iww(1, 2, 0), iww(2, 2, 0)};
}

template <>
void split_face_indices<3, Side::FRONT>(AmrFaces& faces) {
    faces[Side::FRONT0].vertices = {iww(0, 0, 2), iww(1, 0, 2), iww(0, 1, 2), iww(1, 1, 2)};
    faces[Side::FRONT1].vertices = {iww(1, 0, 2), iww(2, 0, 2), iww(1, 1, 2), iww(2, 1, 2)};
    faces[Side::FRONT2].vertices = {iww(0, 1, 2), iww(1, 1, 2), iww(0, 2, 2), iww(1, 2, 2)};
    faces[Side::FRONT3].vertices = {iww(1, 1, 2), iww(2, 1, 2), iww(1, 2, 2), iww(2, 2, 2)};
}

/// @brief Устанавливает нормали и площади новых подграней
template <int dim>
void setup_face_features(AmrFace& face, AmrVertices& vertices);

template <>
void setup_face_features<2>(AmrFace& face, AmrVertices& vertices) {
    // Точка внутри ячейки
    Vector3d C = vertices[iww(1, 1)];

    // Вершины грани
    ShortList1D vl = {
            (Vector3d &) vertices[face.vertices[0]],
            (Vector3d &) vertices[face.vertices[1]],
    };

    // Установить длину и нормаль
    face.area = geom::length(vl);
    face.normal = geom::normal(vl, C);
}

template <>
void setup_face_features<3>(AmrFace& face, AmrVertices& vertices) {
    // Точка внутри ячейки
    Vector3d C = (Vector3d&)vertices[iww(1, 1, 1)];

    // Вершины подграней
    ShortList2D vl = {
            (Vector3d &) vertices[face.vertices[0]],
            (Vector3d &) vertices[face.vertices[1]],
            (Vector3d &) vertices[face.vertices[2]],
            (Vector3d &) vertices[face.vertices[3]],
    };

    // Установить площадь и нормаль
    face.area = geom::area(vl);
    face.normal = geom::normal(vl, C);
}

/// @brief Устанавливает нормали и площади новых подграней
template <int dim>
void split_face_features(AmrFaces& faces, AmrVertices& vertices, int side) {
    setup_face_features<dim>(faces[side], vertices);
    setup_face_features<dim>(faces[side + 6], vertices);

    if (dim > 2) {
        setup_face_features<dim>(faces[side + 12], vertices);
        setup_face_features<dim>(faces[side + 18], vertices);
    }
}

/// @brief Разбивает грань на стороне side на подграни.
/// Устанавливает площади и нормали новых граней.
/// Размерность и сторона являются аргументами шаблона
template <int dim, int side>
void split_face(AmrFaces& faces, AmrVertices& vertices) {
#if SCRUTINY
    if (faces[side + 6].is_actual()) {
        throw std::runtime_error("Try to cut complex face");
    }
#endif

    split_face_prepare<dim, side>(faces);

    split_face_indices<dim, side>(faces);

    split_face_features<dim>(faces, vertices, side);
}

/// @brief Разбивает грань на стороне side на подграни.
/// Устанавливает площади и нормали новых граней.
/// Преобразует аргумент функции side в аргумент шаблона.
template <int dim>
void split_face(AmrFaces& faces, AmrVertices& vertices, int side);

template <>
void split_face<2>(AmrFaces& faces, AmrVertices& vertices, int side) {
    switch (side) {
        case Side::LEFT:
            split_face<2, Side::L>(faces, vertices);
            break;
        case Side::RIGHT:
            split_face<2, Side::R>(faces, vertices);
            break;
        case Side::BOTTOM:
            split_face<2, Side::B>(faces, vertices);
            break;
        default:
            split_face<2, Side::T>(faces, vertices);
            break;
    }
}

template <>
void split_face<3>(AmrFaces& faces, AmrVertices& vertices, int side) {
    switch (side) {
        case Side::LEFT:
            split_face<3, Side::L>(faces, vertices);
            break;
        case Side::RIGHT:
            split_face<3, Side::R>(faces, vertices);
            break;
        case Side::BOTTOM:
            split_face<3, Side::B>(faces, vertices);
            break;
        case Side::TOP:
            split_face<3, Side::T>(faces, vertices);
            break;
        case Side::BACK:
            split_face<3, Side::X>(faces, vertices);
            break;
        default:
            split_face<3, Side::F>(faces, vertices);
            break;
    }
}

/// @brief Разбивает грань на стороне side на подграни.
/// Устанавливает площади и нормали новых граней.
template <int dim>
void split_face(AmrCell& cell, int side) {
    split_face<dim>(cell.faces, cell.vertices, side);
}

/// @brief Объединяет подграни в одну грань на стороне side
template<int dim, Side side>
void merge_faces(AmrCell& cell) {
    auto &faces = cell.faces;
    auto &vertices = cell.vertices;

    faces[side].vertices = face_indices<dim, side>();

    faces[side + 6].set_undefined();

    if (dim > 2) {
        faces[side + 12].set_undefined();
        faces[side + 18].set_undefined();
    }

    setup_face_features<dim>(faces[side], vertices);
}

/// @brief Объединяет подграни в одну грань на стороне side
template<int dim>
void merge_faces(AmrCell& cell, Side side);

template<>
void merge_faces<2>(AmrCell& cell, Side side) {
    switch (side) {
        case Side::LEFT:
            merge_faces<2, Side::L>(cell);
            break;
        case Side::RIGHT:
            merge_faces<2, Side::R>(cell);
            break;
        case Side::BOTTOM:
            merge_faces<2, Side::B>(cell);
            break;
        default:
            merge_faces<2, Side::T>(cell);
            break;
    }
}

template<>
void merge_faces<3>(AmrCell& cell, Side side) {
    switch (side) {
        case Side::LEFT:
            merge_faces<3, Side::L>(cell);
            break;
        case Side::RIGHT:
            merge_faces<3, Side::R>(cell);
            break;
        case Side::BOTTOM:
            merge_faces<3, Side::B>(cell);
            break;
        case Side::TOP:
            merge_faces<3, Side::T>(cell);
            break;
        case Side::BACK:
            merge_faces<3, Side::X>(cell);
            break;
        default:
            merge_faces<3, Side::F>(cell);
            break;
    }
}

} // namespace amr
} // namespace mesh
} // namespace zephyr