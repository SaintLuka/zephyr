/// @brief Файл содержит операции разбиения и огрубления для граней ячейки.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#include <zephyr/geom/primitives/line.h>
#include <zephyr/mesh/primitives/side.h>
#include <zephyr/mesh/euler/amr_cells.h>
#include <zephyr/mesh/amr/common.h>


namespace zephyr::mesh::amr {

/// @brief Устанавливает на новых подгранях значение adjacent и
/// boundary со старой грани
template <int dim, int int_side>
void split_face_prepare(index_t iface, AmrFaces& faces) {
    constexpr Side<dim> side{int_side};
    
    faces.boundary[iface + side[1]] = faces.boundary[iface + side];
    faces.adjacent.rank [iface + side[1]] = faces.adjacent.rank [iface + side];
    faces.adjacent.index[iface + side[1]] = faces.adjacent.index[iface + side];
    faces.adjacent.alien[iface + side[1]] = faces.adjacent.alien[iface + side];
    faces.adjacent.basic[iface + side[1]] = faces.adjacent.basic[iface + side];

    if constexpr (dim > 2) {
        faces.boundary[iface + side[2]] = faces.boundary[iface + side];
        faces.adjacent.rank [iface + side[2]] = faces.adjacent.rank [iface + side];
        faces.adjacent.index[iface + side[2]] = faces.adjacent.index[iface + side];
        faces.adjacent.alien[iface + side[2]] = faces.adjacent.alien[iface + side];
        faces.adjacent.basic[iface + side[2]] = faces.adjacent.basic[iface + side];

        faces.boundary[iface + side[3]] = faces.boundary[iface + side];
        faces.adjacent.rank [iface + side[3]] = faces.adjacent.rank [iface + side];
        faces.adjacent.index[iface + side[3]] = faces.adjacent.index[iface + side];
        faces.adjacent.alien[iface + side[3]] = faces.adjacent.alien[iface + side];
        faces.adjacent.basic[iface + side[3]] = faces.adjacent.basic[iface + side];
    }
}

/// @brief Устанавливает на новых подгранях индексы вершин.
template <int dim, int int_side>
void split_face_indices(index_t iface, AmrFaces& faces) {
    constexpr Side<dim> side{int_side};

    faces.vertices[iface + side] = side.cf();
    faces.vertices[iface + side[1]] = side[1].cf();

    if constexpr (dim > 2) {
        faces.vertices[iface + side[2]] = side[2].cf();
        faces.vertices[iface + side[3]] = side[3].cf();
    }
}

/// @brief Устанавливает нормали и площади новых подграней
template <int dim>
void setup_face_features(index_t iface, AmrFaces& faces, const SqMap<dim>& vertices, bool axial = false);

template <>
inline void setup_face_features<2>(index_t iface, AmrFaces& faces, const SqQuad& vertices, bool axial) {
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
inline void setup_face_features<3>(index_t iface, AmrFaces& faces, const SqCube& vertices, bool axial) {
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
template <int dim, int side_in>
void split_face_features(index_t iface, AmrFaces& faces,
                         const SqMap<dim>& vertices, bool axial) {
    constexpr Side<dim> side{side_in};

    setup_face_features<dim>(iface + side, faces, vertices, axial);
    setup_face_features<dim>(iface + side[1], faces, vertices, axial);

    if constexpr (dim > 2) {
        setup_face_features<dim>(iface + side[2], faces, vertices);
        setup_face_features<dim>(iface + side[3], faces, vertices);
    }
}

/// @brief Разбивает грань на стороне side на подграни.
/// Устанавливает площади и нормали новых граней.
/// Размерность и сторона являются аргументами шаблона
template <int dim, int side>
void split_face(index_t iface, AmrFaces& faces, const SqMap<dim>& vertices, bool axial = false) {
#if SCRUTINY
    if (faces.is_actual(iface + Side<dim>(side)[1])) {
        throw std::runtime_error("Try to cut complex face");
    }
#endif

    split_face_prepare<dim, side>(iface, faces);

    split_face_indices<dim, side>(iface, faces);

    split_face_features<dim, side>(iface, faces, vertices, axial);
}

/// @brief Разбивает грань на стороне side на подграни.
/// Устанавливает площади и нормали новых граней.
/// Преобразует аргумент функции side в аргумент шаблона.
template <int dim>
void split_face(index_t iface, AmrFaces& faces, const SqMap<dim>& vertices, Side<dim> side, bool axial);

template <>
inline void split_face<2>(index_t iface, AmrFaces& faces, const SqQuad& vertices, Side2D side, bool axial) {
    switch (side) {
        case Side2D::LEFT:
            split_face<2, Side2D::L>(iface, faces, vertices, axial);
            break;
        case Side2D::RIGHT:
            split_face<2, Side2D::R>(iface, faces, vertices, axial);
            break;
        case Side2D::BOTTOM:
            split_face<2, Side2D::B>(iface, faces, vertices, axial);
            break;
        default:
            split_face<2, Side2D::T>(iface, faces, vertices, axial);
            break;
    }
}

template <>
inline void split_face<3>(index_t iface, AmrFaces& faces, const SqCube& vertices, Side3D side, bool axial) {
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
void split_face(index_t ic, AmrCells& cells, int side) {
    split_face<dim>(cells.face_begin[ic], cells.faces, cells.mapping<dim>(ic), side, cells.axial());
}

/// @brief Объединяет подграни в одну грань на стороне side
template<int dim, int int_side>
void merge_faces(index_t ic, AmrCells& cells) {
    constexpr Side<dim> side{int_side};

    const index_t iface  = cells.face_begin[ic];
    const SqMap<dim> &vertices = cells.mapping<dim>(ic);

    cells.faces.vertices[iface + side] = side.sf();

    cells.faces.set_undefined(iface + side[1]);

    if constexpr (dim > 2) {
        cells.faces.set_undefined(iface + side[2]);
        cells.faces.set_undefined(iface + side[3]);
    }

    setup_face_features<dim>(iface + side, cells.faces, vertices, cells.axial());
}

/// @brief Объединяет подграни в одну грань на стороне side
template<int dim>
void merge_faces(index_t ic, AmrCells& cells, Side<dim> side);

template<>
inline void merge_faces<2>(index_t ic, AmrCells& cells, Side2D side) {
    switch (side) {
        case Side2D::LEFT:
            merge_faces<2, Side2D::L>(ic, cells);
            break;
        case Side2D::RIGHT:
            merge_faces<2, Side2D::R>(ic, cells);
            break;
        case Side2D::BOTTOM:
            merge_faces<2, Side2D::B>(ic, cells);
            break;
        default:
            merge_faces<2, Side2D::T>(ic, cells);
            break;
    }
}

template<>
inline void merge_faces<3>(index_t ic, AmrCells& cells, Side3D side) {
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

} // namespace zephyr::mesh::amr