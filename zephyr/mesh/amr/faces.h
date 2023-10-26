/// @file Файл содержит операции разбиения и огрубления для граней ячейки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/geom/primitives/side.h>
#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/mesh/amr/common.h>

namespace zephyr::mesh::amr {

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
void split_face_indices(AmrFaces& faces) {
    faces[side].vertices = face_indices::cf<dim, Side(side)>();
    faces[side + 6].vertices = face_indices::cf<dim, Side(side + 6)>();

    if (dim > 2) {
        faces[side + 12].vertices = face_indices::cf<dim, Side(side + 12)>();
        faces[side + 18].vertices = face_indices::cf<dim, Side(side + 18)>();
    }
}

/// @brief Устанавливает нормали и площади новых подграней
template <int dim>
void setup_face_features(AmrFace& face, SqCube& vertices);

template <>
void setup_face_features<2>(AmrFace& face, SqCube& vertices) {
    // Точка внутри ячейки
    Vector3d C = vertices.vs<0, 0>();

    // Вершины грани
    Line vl = {
            vertices[face.vertices[0]],
            vertices[face.vertices[1]],
    };

    // Установить длину и нормаль
    face.area   = vl.length();
    face.center = vl.center();
    face.normal = vl.normal(C);
}

template <>
void setup_face_features<3>(AmrFace& face, SqCube& vertices) {
    // Точка внутри ячейки
    const Vector3d& C = vertices.vs<0, 0, 0>();

    // Вершины подграней
    Quad vl = {
            vertices[face.vertices[0]],
            vertices[face.vertices[1]],
            vertices[face.vertices[2]],
            vertices[face.vertices[3]],
    };

    // Установить площадь и нормаль
    face.area   = vl.area();
    face.center = vl.center();
    face.normal = vl.normal(C);
}

/// @brief Устанавливает нормали и площади новых подграней
template <int dim>
void split_face_features(AmrFaces& faces, SqCube& vertices, int side) {
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
void split_face(AmrFaces& faces, SqCube& vertices) {
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
void split_face(AmrFaces& faces, SqCube& vertices, int side);

template <>
void split_face<2>(AmrFaces& faces, SqCube& vertices, int side) {
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
void split_face<3>(AmrFaces& faces, SqCube& vertices, int side) {
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

    faces[side].vertices = face_indices::sf<dim, side>();

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

} // namespace zephyr::mesh::amr