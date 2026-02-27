#pragma once
#include <string>

#include <zephyr/mesh/side.h>

namespace zephyr::geom::indexing {

constexpr int triangle_face_size[3] = {2, 2, 2};
constexpr int triangle_face_indices[3][2] = {
    {0, 1}, {1, 2}, {2, 0}
};

constexpr int quad_face_size[4] = {2, 2, 2, 2};
constexpr int quad_face_indices[4][2] = {
    {0, 3}, {1, 2}, {0, 1}, {3, 2}
};

constexpr int tetra_face_size[4] = {3, 3, 3, 3};
constexpr int tetra_face_indices[5][4] = {
    {0, 2, 1},
    {0, 1, 3},
    {1, 2, 3},
    {0, 3, 2}
};

constexpr int pyramid_face_size[4] = {4, 3, 3, 3};
constexpr int pyramid_face_indices[5][4] = {
    {3, 2, 1, 0},
    {0, 1, 4},
    {1, 2, 4},
    {2, 3, 4},
    {0, 4, 3}
};

constexpr int wedge_face_size[5] = {4, 4, 4, 3, 3};
constexpr int wedge_face_indices[5][4] =
    {{0, 3, 4, 1},
     {1, 4, 5, 2},
     {0, 2, 5, 3},
     {0, 1, 2},
     {3, 5, 4}};

constexpr int hexahedron_face_size[6] = {4, 4, 4, 4, 4, 4};
constexpr int hexahedron_face_indices[6][4] = {
    {0, 3, 2, 1},
    {4, 5, 6, 7},
    {0, 3, 4, 7},
    {1, 2, 6, 5},
    {0, 1, 5, 4},
    {2, 3, 7, 6}
};

// Индексация вершин в кубической AMR-ячейке. Повторяет индексацию функции SqCube::iss
template<int i, int j, int k = -1>
static constexpr int iss() {
    static_assert(i * i <= 1 && j * j <= 1 && k * k <= 1, "Available indices: {-1, 0, 1}");
    return 9 * (k + 1) + 3 * (j + 1) + i + 1;
}

constexpr int amr2d_face_size[4] = {2, 2, 2, 2};
constexpr int amr2d_face_indices[4][2] = {
    {iss<-1, -1>(), iss<-1, +1>()},
    {iss<+1, -1>(), iss<+1, +1>()},
    {iss<-1, -1>(), iss<+1, -1>()},
    {iss<-1, +1>(), iss<+1, +1>()},
};

constexpr int amr3d_face_size[6] = {4, 4, 4, 4};
constexpr int amr3d_face_indices[6][4] = {
    {iss<-1, -1, -1>(), iss<-1, +1, -1>(), iss<-1, -1, +1>(), iss<-1, +1, +1>()},
    {iss<+1, -1, -1>(), iss<+1, +1, -1>(), iss<+1, -1, +1>(), iss<+1, +1, +1>()},
    {iss<-1, -1, -1>(), iss<+1, -1, -1>(), iss<-1, -1, +1>(), iss<+1, -1, +1>()},
    {iss<-1, +1, -1>(), iss<+1, +1, -1>(), iss<-1, +1, +1>(), iss<+1, +1, +1>()},
    {iss<-1, -1, -1>(), iss<+1, -1, -1>(), iss<-1, +1, -1>(), iss<+1, +1, -1>()},
    {iss<-1, -1, +1>(), iss<+1, -1, +1>(), iss<-1, +1, +1>(), iss<+1, +1, +1>()},
};

/// @brief Адаптивный элемент?
constexpr bool is_adaptive(CellType type) {
    return type == CellType::AMR2D || type == CellType::AMR3D;
}

/// @brief Динамическое число узлов и граней?
constexpr bool is_dynamic(CellType type) {
    return type == CellType::POLYGON || type == CellType::POLYHEDRON;
}

/// @brief Размерность примитива
constexpr int get_dimension(CellType type) {
    if (type == CellType::TRIANGLE ||
        type == CellType::QUAD ||
        type == CellType::POLYGON ||
        type == CellType::AMR2D) {
        return 2;
    }
    if (type == CellType::TETRA  ||
        type == CellType::PYRAMID ||
        type == CellType::WEDGE ||
        type == CellType::HEXAHEDRON ||
        type == CellType::POLYHEDRON ||
        type == CellType::AMR3D) {
        return 3;
    }
    throw std::runtime_error("dimension(CellType) error: unknown cell type");
}

/// @brief Число узлов элемента
constexpr int n_nodes(CellType type) {
    switch (type) {
        case CellType::TRIANGLE:   return 3;
        case CellType::QUAD:       return 4;
        case CellType::AMR2D:      return 9;
        case CellType::TETRA:      return 4;
        case CellType::PYRAMID:    return 5;
        case CellType::WEDGE:      return 6;
        case CellType::HEXAHEDRON: return 8;
        case CellType::AMR3D:      return 27;
        default: return -1; // dynamic
    }
}

/// @brief Число граней элемента
constexpr int n_faces(CellType type) {
    switch (type) {
        case CellType::TRIANGLE:   return 3;
        case CellType::QUAD:
        case CellType::AMR2D:
        case CellType::TETRA:      return 4;
        case CellType::PYRAMID:
        case CellType::WEDGE:      return 5;
        case CellType::HEXAHEDRON:
        case CellType::AMR3D:      return 6;
        default: return -1; // dynamic
    }
}

/// @brief Число узлов грани
constexpr int face_size(CellType type, int side) {
    switch (type) {
        case CellType::TRIANGLE:
        case CellType::QUAD:
        case CellType::POLYGON:
        case CellType::AMR2D:      return 2;
        case CellType::TETRA:      return 3;
        case CellType::PYRAMID:    return pyramid_face_size[side];
        case CellType::WEDGE:      return wedge_face_size[side];
        case CellType::HEXAHEDRON:
        case CellType::AMR3D:      return 4;
        default: return -1; // dynamic
    }
}

/// @brief Число узлов грани
constexpr int face_size(CellType type, int side, int idx) {
    z_assert(idx < face_size(type, side), "bad face index");
    switch (type) {
        case CellType::TRIANGLE:   return triangle_face_indices[side][idx];
        case CellType::QUAD:       return quad_face_indices[side][idx];
        case CellType::AMR2D:      return amr2d_face_indices[side][idx];
        case CellType::TETRA:      return tetra_face_indices[side][idx];
        case CellType::PYRAMID:    return pyramid_face_indices[side][idx];
        case CellType::WEDGE:      return wedge_face_indices[side][idx];
        case CellType::HEXAHEDRON: return hexahedron_face_indices[side][idx];
        case CellType::AMR3D:      return amr3d_face_indices[side][idx];
        default: return -1; // dynamic
    }
}

constexpr std::string to_string(CellType type) {
    switch (type) {
        case CellType::TRIANGLE:   return "Triangle";
        case CellType::QUAD:       return "Quad";
        case CellType::POLYGON:    return "Polygon";
        case CellType::AMR2D:      return "AMR2D";
        case CellType::TETRA:      return "Tetrahedron";
        case CellType::PYRAMID:    return "Pyramid";
        case CellType::WEDGE:      return "Wedge";
        case CellType::HEXAHEDRON: return "Hexahedron";
        case CellType::AMR3D:      return "AMR3D";
        case CellType::POLYHEDRON: return "Polyhedron";
        default: return "Unknown";
    }
}

} // namespace zephyr::geom::indexing