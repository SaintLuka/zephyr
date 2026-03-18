#pragma once
#include <string>
#include <algorithm>

#include <zephyr/configuration.h>
#include <zephyr/geom/side.h>
#include <zephyr/geom/cell_type.h>

namespace zephyr::geom::indexing {

/// @brief Число вершин простой грани.
constexpr int VpF(int dim) { return dim < 3 ? 2 : 4; }

/// @brief Количество подграней.
constexpr int FpF(int dim) { return dim < 3 ? 2 : 4; }

/// @brief Число граней у простой ячейки.
constexpr int FpC(int dim) { return dim < 3 ? 4 : 6; }

/// @brief Количество дочерних ячеек.
constexpr int CpC(int dim) { return dim < 3 ? 4 : 8; }

// Как хранится в AmrCells
using node_array = std::array<short, 8>;

template <int size, typename container>
constexpr node_array convert(const container& arr) {
    static_assert(size < 8);
    node_array res = {-1, -1, -1, -1, -1, -1, -1, -1};
    std::ranges::copy_n(arr, size, res.begin());
    return res;
}

// --------------------------------- Элементы ---------------------------------

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

/// @brief Тип элемента в строку
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

// ----------------------------------- Узлы -----------------------------------

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

// ----------------------------------------------------------------------------

/// @brief Индексация узлов и граней в треугольнике (лол)
namespace tri {

constexpr int n_faces = 3;
constexpr int face_size[3] = {2, 2, 2};
constexpr int face_nodes[3][2] = {
    {0, 1}, {1, 2}, {2, 0}
};

/// @brief Индексы вершин грани
constexpr node_array sf(int side) { return convert<2>(face_nodes[side]); }

} // namespace indexing::tri

// ----------------------------------------------------------------------------

/// @brief Индексация узлов и граней в простой QUAD-ячейке
namespace quad {

/// @brief Нумерация узлов в простом четырёхугольнике.
/// Против часовой стрелки от 0 до 3 с левого нижнего угла.
template <int i, int j>
constexpr int vs() {
    static_assert(0 <= std::min(i, j) && std::max(i, j) <= 1, "Available indices: {0, 1}");
    return i == 0 ? (j == 0 ? 0 : 3) : (j == 0 ? 1 : 2);
}

constexpr int n_faces = 4;
constexpr int face_size[4] = {2, 2, 2, 2};
constexpr int face_nodes[4][2] = {
    {vs<0,0>(), vs<0,1>()}, // L
    {vs<1,0>(), vs<1,1>()}, // R
    {vs<0,0>(), vs<1,0>()}, // B
    {vs<0,1>(), vs<1,1>()}, // T
};

/// @brief Индексы вершин грани
constexpr node_array sf(int side) { return convert<2>(face_nodes[side]); }

} // namespace indexing::quad

// ----------------------------------------------------------------------------

namespace poly {

constexpr node_array sf(int count, int side) {
    node_array vertices;
    std::ranges::fill(vertices, -1);
    vertices[0] = static_cast<short>(side);
    vertices[1] = static_cast<short>((side + 1) % count);
    return vertices;
}

} // namespace indexing::poly

// ----------------------------------------------------------------------------

namespace tetra {

constexpr int n_nodes = 4;
constexpr int n_faces = 4;
constexpr int face_size[4] = {3, 3, 3, 3};
constexpr int face_nodes[4][3] = {
    {0, 2, 1},
    {0, 1, 3},
    {1, 2, 3},
    {0, 3, 2}
};

/// @brief Индексы вершин грани
constexpr node_array sf(int side) { return convert<3>(face_nodes[side]); }

} // namespace indexing::tetra

// ----------------------------------------------------------------------------

namespace pyramid {

constexpr int n_nodes = 5;
constexpr int n_faces = 5;
constexpr int face_size[5] = {4, 3, 3, 3, 3};
constexpr int face_nodes[5][4] = {
    {3, 2, 1, 0},
    {0, 1, 4},
    {1, 2, 4},
    {2, 3, 4},
    {0, 4, 3}
};

/// @brief Индексы вершин грани
constexpr node_array sf(int side) {
    if (face_size[side] == 3) {
        return convert<3>(face_nodes[side]);
    }
    else {
        return convert<4>(face_nodes[side]);
    }
}

} // namespace indexing::pyramid

// ----------------------------------------------------------------------------

namespace wedge {

/// @brief Нумерация узлов в треугольной призме.
/// b - индекс основания, i - индекс вершины основания
template <int b, int i>
constexpr int vs() {
    static_assert(0 <= b && b <= 1 && 0 <= i && i <= 2, "Available indices: {0, 1, 2}");
    return 3 * b + i;
}

constexpr int n_nodes = 6;
constexpr int n_faces = 5;
constexpr int face_size[5] = {4, 4, 4, 3, 3};
constexpr int face_nodes[5][4] ={
    {0, 3, 4, 1},
    {1, 4, 5, 2},
    {0, 2, 5, 3},
    {0, 1, 2},
    {3, 5, 4}
};

/// @brief Индексы вершин грани
constexpr node_array sf(int side) {
    if (face_size[side] == 3) {
        return convert<3>(face_nodes[side]);
    }
    else {
        return convert<4>(face_nodes[side]);
    }
}

} // namespace indexing::wedge

// ----------------------------------------------------------------------------

namespace hex {

/// @brief Нумерация узлов в шестиграннике (соглашение VTK).
template <int i, int j, int k>
constexpr int vs() {
    static_assert(0 <= std::min({i, j, k}) && std::max({i, j, k}) <= 1, "Available indices: {0, 1}");
    return quad::vs<i, j>() + 4 * k;
}

constexpr int n_nodes = 8;
constexpr int n_faces = 6;
constexpr int face_size[6] = {4, 4, 4, 4, 4, 4};
constexpr int face_nodes[6][4] = {
    {vs<0,0,0>(), vs<0,0,1>(), vs<0,1,1>(), vs<0,1,0>()}, // L
    {vs<1,0,0>(), vs<1,1,0>(), vs<1,1,1>(), vs<1,0,1>()}, // R
    {vs<0,0,0>(), vs<1,0,0>(), vs<1,0,1>(), vs<0,0,1>()}, // B
    {vs<1,1,0>(), vs<0,1,0>(), vs<0,1,1>(), vs<1,1,1>()}, // T
    {vs<0,0,0>(), vs<0,1,0>(), vs<1,1,0>(), vs<1,0,0>()}, // X
    {vs<0,0,1>(), vs<1,0,1>(), vs<1,1,1>(), vs<0,1,1>()}, // F
};

/// @brief Индексы вершин грани
constexpr node_array sf(int side) { return convert<4>(face_nodes[side]); }

} // namespace indexing::hex

// ----------------------------------------------------------------------------

/// @brief Индексация узлов и граней в AMR-ячейках
namespace amr {

/// @brief Индексация вершин в кубической AMR-ячейке. Повторяет индексацию функции SqCube::iss
template <int i, int j, int k = -1>
constexpr int vs() {
    static_assert(i * i <= 1 && j * j <= 1 && k * k <= 1, "Available indices: {-1, 0, 1}");
    return 9 * (k + 1) + 3 * (j + 1) + i + 1;
}

constexpr int n_faces2d = 4;
constexpr int face_size2d[4] = {2, 2, 2, 2};
constexpr int face_nodes2d[4][2] = {
    {vs<-1, -1>(), vs<-1, +1>()},
    {vs<+1, -1>(), vs<+1, +1>()},
    {vs<-1, -1>(), vs<+1, -1>()},
    {vs<-1, +1>(), vs<+1, +1>()},
};

constexpr int sqface_size2d[4] = {3, 3, 3, 3};
constexpr int sqface_nodes2d[4][3] = {
    {vs<-1, -1>(), vs<-1, 0>(), vs<-1, +1>()},
    {vs<+1, -1>(), vs<+1, 0>(), vs<+1, +1>()},
    {vs<-1, -1>(), vs<0, -1>(), vs<+1, -1>()},
    {vs<-1, +1>(), vs<0, +1>(), vs<+1, +1>()},
};

constexpr int n_faces3d = 6;
constexpr int face_size3d[6] = {4, 4, 4, 4};
constexpr int face_nodes3d[6][4] = {
    {vs<-1, -1, -1>(), vs<-1, +1, -1>(), vs<-1, -1, +1>(), vs<-1, +1, +1>()},
    {vs<+1, -1, -1>(), vs<+1, +1, -1>(), vs<+1, -1, +1>(), vs<+1, +1, +1>()},
    {vs<-1, -1, -1>(), vs<+1, -1, -1>(), vs<-1, -1, +1>(), vs<+1, -1, +1>()},
    {vs<-1, +1, -1>(), vs<+1, +1, -1>(), vs<-1, +1, +1>(), vs<+1, +1, +1>()},
    {vs<-1, -1, -1>(), vs<+1, -1, -1>(), vs<-1, +1, -1>(), vs<+1, +1, -1>()},
    {vs<-1, -1, +1>(), vs<+1, -1, +1>(), vs<-1, +1, +1>(), vs<+1, +1, +1>()},
};

constexpr int sqface_size3d[6] = {9, 9, 9, 9};
constexpr int sqface_nodes3d[6][9] = {
    {vs<-1, -1, -1>(), vs<-1,  0, -1>(), vs<-1, +1, -1>(),
     vs<-1, -1,  0>(), vs<-1,  0,  0>(), vs<-1, +1,  0>(),
     vs<-1, -1, +1>(), vs<-1,  0, +1>(), vs<-1, +1, +1>()}, // LEFT
    {vs<+1, -1, -1>(), vs<+1,  0, -1>(), vs<+1, +1, -1>(),
     vs<+1, -1,  0>(), vs<+1,  0,  0>(), vs<+1, +1,  0>(),
     vs<+1, -1, +1>(), vs<+1,  0, +1>(), vs<+1, +1, +1>()}, // RIGHT
    {vs<-1, -1, -1>(), vs< 0, -1, -1>(), vs<+1, -1, -1>(),
     vs<-1, -1,  0>(), vs< 0, -1,  0>(), vs<+1, -1,  0>(),
     vs<-1, -1, +1>(), vs< 0, -1, +1>(), vs<+1, -1, +1>()}, // BOTTOM
    {vs<-1, +1, -1>(), vs< 0, +1, -1>(), vs<+1, +1, -1>(),
     vs<-1, +1,  0>(), vs< 0, +1,  0>(), vs<+1, +1,  0>(),
     vs<-1, +1, +1>(), vs< 0, +1, +1>(), vs<+1, +1, +1>()}, // TOP
    {vs<-1, -1, -1>(), vs< 0, -1, -1>(), vs<+1, -1, -1>(),
     vs<-1,  0, -1>(), vs< 0,  0, -1>(), vs<+1,  0, -1>(),
     vs<-1, +1, -1>(), vs< 0, +1, -1>(), vs<+1, +1, -1>()}, // BACK
    {vs<-1, -1, +1>(), vs< 0, -1, +1>(), vs<+1, -1, +1>(),
     vs<-1,  0, +1>(), vs< 0,  0, +1>(), vs<+1,  0, +1>(),
     vs<-1, +1, +1>(), vs< 0, +1, +1>(), vs<+1, +1, +1>()}, // FRONT
};

constexpr int n_subfaces2d = 8;
constexpr int subface_nodes2d[8][2] = {
    {vs<-1,-1>(), vs<-1, 0>()}, // L[0]
    {vs<+1,-1>(), vs<+1, 0>()}, // R[0]
    {vs<-1,-1>(), vs< 0,-1>()}, // B[0]
    {vs<-1,+1>(), vs< 0,+1>()}, // T[0]

    {vs<-1, 0>(), vs<-1,+1>()}, // L[1]
    {vs<+1, 0>(), vs<+1,+1>()}, // R[1]
    {vs< 0,-1>(), vs<+1,-1>()}, // B[1]
    {vs< 0,+1>(), vs<+1,+1>()}, // T[1]
};

constexpr int n_subfaces3d = 24;
constexpr int subface_nodes3d[24][4] = {
    {vs<-1,-1,-1>(), vs<-1, 0,-1>(), vs<-1,-1, 0>(), vs<-1, 0, 0>()}, // L[0]
    {vs<+1,-1,-1>(), vs<+1, 0,-1>(), vs<+1,-1, 0>(), vs<+1, 0, 0>()}, // R[0]
    {vs<-1,-1,-1>(), vs< 0,-1,-1>(), vs<-1,-1, 0>(), vs< 0,-1, 0>()}, // B[0]
    {vs<-1,+1,-1>(), vs< 0,+1,-1>(), vs<-1,+1, 0>(), vs< 0,+1, 0>()}, // T[0]
    {vs<-1,-1,-1>(), vs<0, -1,-1>(), vs<-1, 0,-1>(), vs< 0, 0,-1>()}, // Z[0]
    {vs<-1,-1,+1>(), vs< 0,-1,+1>(), vs<-1, 0,+1>(), vs< 0, 0,+1>()}, // F[0]

    {vs<-1, 0,-1>(), vs<-1,+1,-1>(), vs<-1, 0, 0>(), vs<-1,+1, 0>()}, // L[1]
    {vs<+1, 0,-1>(), vs<+1,+1,-1>(), vs<+1, 0, 0>(), vs<+1,+1, 0>()}, // R[1]
    {vs< 0,-1,-1>(), vs<+1,-1,-1>(), vs< 0,-1, 0>(), vs<+1,-1, 0>()}, // B[1]
    {vs< 0,+1,-1>(), vs<+1,+1,-1>(), vs< 0,+1, 0>(), vs<+1,+1, 0>()}, // T[1]
    {vs< 0,-1,-1>(), vs<+1,-1,-1>(), vs< 0, 0,-1>(), vs<+1, 0,-1>()}, // Z[1]
    {vs< 0,-1,+1>(), vs<+1,-1,+1>(), vs< 0, 0,+1>(), vs<+1, 0,+1>()}, // F[1]

    {vs<-1,-1, 0>(), vs<-1, 0, 0>(), vs<-1,-1,+1>(), vs<-1, 0,+1>()}, // L[2]
    {vs<+1,-1, 0>(), vs<+1, 0, 0>(), vs<+1,-1,+1>(), vs<+1, 0,+1>()}, // R[2]
    {vs<-1,-1, 0>(), vs< 0,-1, 0>(), vs<-1,-1,+1>(), vs< 0,-1,+1>()}, // B[2]
    {vs<-1,+1, 0>(), vs< 0,+1, 0>(), vs<-1,+1,+1>(), vs< 0,+1,+1>()}, // T[2]
    {vs<-1, 0,-1>(), vs< 0, 0,-1>(), vs<-1,+1,-1>(), vs< 0,+1,-1>()}, // Z[2]
    {vs<-1, 0,+1>(), vs< 0, 0,+1>(), vs<-1,+1,+1>(), vs< 0,+1,+1>()}, // F[2]

    {vs<-1, 0, 0>(), vs<-1,+1, 0>(), vs<-1, 0,+1>(), vs<-1,+1,+1>()}, // L[3]
    {vs<+1, 0, 0>(), vs<+1,+1, 0>(), vs<+1, 0,+1>(), vs<+1,+1,+1>()}, // R[3]
    {vs< 0,-1, 0>(), vs<+1,-1, 0>(), vs< 0,-1,+1>(), vs<+1,-1,+1>()}, // B[3]
    {vs< 0,+1, 0>(), vs<+1,+1, 0>(), vs< 0,+1,+1>(), vs<+1,+1,+1>()}, // T[3]
    {vs< 0, 0,-1>(), vs<+1, 0,-1>(), vs< 0,+1,-1>(), vs<+1,+1,-1>()}, // Z[3]
    {vs< 0, 0,+1>(), vs<+1, 0,+1>(), vs< 0,+1,+1>(), vs<+1,+1,+1>()}, // F[3]
};

/// @brief Индексы вершин простой грани в AMR-ячейке (sf -- simple face)
constexpr node_array sf(Side2D side) { return convert<2>(face_nodes2d[side]); }

/// @brief Индексы вершин простой грани в AMR-ячейке (sf -- simple face)
constexpr node_array sf(Side3D side) { return convert<4>(face_nodes3d[side]); }

/// @brief Индексы вершин части составной грани AMR-ячейки (cf -- complex face)
constexpr node_array cf(Side2D side) { return convert<2>(subface_nodes2d[side]); }

/// @brief Индексы вершин части составной грани AMR-ячейки (cf -- complex face)
constexpr node_array cf(Side3D side) { return convert<4>(subface_nodes3d[side]); }

} // namespace indexing::amr

// ---------------------------------- Грани -----------------------------------



/// @brief Число граней элемента
constexpr int n_faces(CellType type) {
    switch (type) {
        case CellType::TRIANGLE:   return tri::n_faces;
        case CellType::QUAD:       return quad::n_faces;
        case CellType::AMR2D:      return amr::n_faces2d;
        case CellType::TETRA:      return tetra::n_faces;
        case CellType::PYRAMID:    return pyramid::n_faces;
        case CellType::WEDGE:      return wedge::n_faces;
        case CellType::HEXAHEDRON: return hex::n_faces;
        case CellType::AMR3D:      return amr::n_faces3d;
        default: return -1; // dynamic
    }
}

/// @brief Число узлов грани
constexpr int face_size(CellType type, int side) {
    switch (type) {
        case CellType::TRIANGLE:
        case CellType::QUAD:
        case CellType::POLYGON:    return 2;
        case CellType::AMR2D:      return amr::face_size2d[side];
        case CellType::TETRA:      return tetra::face_size[side];
        case CellType::PYRAMID:    return pyramid::face_size[side];
        case CellType::WEDGE:      return wedge::face_size[side];
        case CellType::HEXAHEDRON: return hex::face_size[side];
        case CellType::AMR3D:      return amr::face_size3d[side];
        default: return -1; // dynamic
    }
}

/// @brief Число узлов грани
constexpr int face_index(CellType type, int side, int idx) {
    z_assert(idx < face_size(type, side), "bad face index");
    switch (type) {
        case CellType::TRIANGLE:   return tri::face_nodes[side][idx];
        case CellType::QUAD:       return quad::face_nodes[side][idx];
        case CellType::AMR2D:      return amr::face_nodes2d[side][idx];
        case CellType::TETRA:      return tetra::face_nodes[side][idx];
        case CellType::PYRAMID:    return pyramid::face_nodes[side][idx];
        case CellType::WEDGE:      return wedge::face_nodes[side][idx];
        case CellType::HEXAHEDRON: return hex::face_nodes[side][idx];
        case CellType::AMR3D:      return amr::face_nodes3d[side][idx];
        default: return -1; // dynamic
    }
}

// ----------------------------- Адаптивная жесть -----------------------------

// К примеру, child_by_face[Side2D::LEFT[1]]
static constexpr short child_by_subface_2D[8] = {
    0, 1, 0, 2, // L[0], R[0], B[0], T[0]
    2, 3, 1, 3, // L[1], R[1], B[1], T[1]
};

// К примеру, children_by_face[Side2D::LEFT]
static constexpr std::array<short, 2> children_by_face_2D[4] = {
    {0, 2},
    {1, 3},
    {0, 1},
    {2, 3},
};

// К примеру, child_by_face[Side3D::LEFT[3]]
static constexpr short child_by_subface_3D[24] = {
    0, 1, 0, 2, 0, 4, // L[0], R[0], B[0], T[0], Z[0], F[0]
    2, 3, 1, 3, 1, 5, // L[1], R[1], B[1], T[1], Z[1], F[1]
    4, 5, 4, 6, 2, 6,
    6, 7, 5, 7, 3, 7,
};

// К примеру, children_by_face[Side3D::LEFT]
static constexpr std::array<short, 4> children_by_face_3D[6] = {
    {0, 2, 4, 6},
    {1, 3, 5, 7},
    {0, 1, 4, 5},
    {2, 3, 6, 7},
    {0, 1, 2, 3},
    {4, 5, 6, 7},
};

/// @brief Получить subface по стороне и индексу дочерней ячейки (частичная генерация)
static constexpr int subface_by_side_child_2D[4][4] = {
    {0, -1, 4, -1},
    {-1, 1, -1, 5},
    {2, 6, -1, -1},
    {-1, -1, 3, 7},
};

/// @brief Получить subface по стороне и индексу дочерней ячейки (частичная генерация)
static constexpr int subface_by_side_child_3D[6][8] = {
    {0, -1, 6, -1, 12, -1, 18, -1},
    {-1, 1, -1, 7, -1, 13, -1, 19},
    {2, 8, -1, -1, 14, 20, -1, -1},
    {-1, -1, 3, 9, -1, -1, 15, 21},
    {4, 10, 16, 22, -1, -1, -1, -1},
    {-1, -1, -1, -1, 5, 11, 17, 23},
};

// [symmetry][subface]
static constexpr int neib_child_by_symm_subface_2D[8][8] = {
    {1, 0, 2, 0, 3, 2, 3, 1},
    {0, 3, 1, 3, 2, 1, 2, 0},
    {3, 2, 0, 2, 1, 0, 1, 3},
    {2, 1, 3, 1, 0, 3, 0, 2},
    {1, 2, 0, 2, 3, 0, 3, 1},
    {2, 3, 1, 3, 0, 1, 0, 2},
    {3, 0, 2, 0, 1, 2, 1, 3},
    {0, 1, 3, 1, 2, 3, 2, 0},
};

// [symmetry][subface]
static constexpr int neib_child_by_symm_subface_3D[48][24] = {
    {1, 0, 2, 0, 4, 0, 3, 2, 3, 1, 5, 1, 5, 4, 6, 4, 6, 2, 7, 6, 7, 5, 7, 3},
    {1, 0, 4, 0, 2, 0, 5, 4, 5, 1, 3, 1, 3, 2, 6, 2, 6, 4, 7, 6, 7, 3, 7, 5},
    {2, 0, 1, 0, 4, 0, 3, 1, 3, 2, 6, 2, 6, 4, 5, 4, 5, 1, 7, 5, 7, 6, 7, 3},
    {4, 0, 1, 0, 2, 0, 5, 1, 5, 4, 6, 4, 6, 2, 3, 2, 3, 1, 7, 3, 7, 6, 7, 5},
    {2, 0, 4, 0, 1, 0, 6, 4, 6, 2, 3, 2, 3, 1, 5, 1, 5, 4, 7, 5, 7, 3, 7, 6},
    {4, 0, 2, 0, 1, 0, 6, 2, 6, 4, 5, 4, 5, 1, 3, 1, 3, 2, 7, 3, 7, 5, 7, 6},
    {0, 1, 3, 1, 5, 1, 2, 3, 2, 0, 4, 0, 4, 5, 7, 5, 7, 3, 6, 7, 6, 4, 6, 2},
    {0, 1, 5, 1, 3, 1, 4, 5, 4, 0, 2, 0, 2, 3, 7, 3, 7, 5, 6, 7, 6, 2, 6, 4},
    {0, 2, 3, 2, 6, 2, 1, 3, 1, 0, 4, 0, 4, 6, 7, 6, 7, 3, 5, 7, 5, 4, 5, 1},
    {0, 4, 5, 4, 6, 4, 1, 5, 1, 0, 2, 0, 2, 6, 7, 6, 7, 5, 3, 7, 3, 2, 3, 1},
    {0, 2, 6, 2, 3, 2, 4, 6, 4, 0, 1, 0, 1, 3, 7, 3, 7, 6, 5, 7, 5, 1, 5, 4},
    {0, 4, 6, 4, 5, 4, 2, 6, 2, 0, 1, 0, 1, 5, 7, 5, 7, 6, 3, 7, 3, 1, 3, 2},
    {3, 1, 0, 1, 5, 1, 2, 0, 2, 3, 7, 3, 7, 5, 4, 5, 4, 0, 6, 4, 6, 7, 6, 2},
    {5, 1, 0, 1, 3, 1, 4, 0, 4, 5, 7, 5, 7, 3, 2, 3, 2, 0, 6, 2, 6, 7, 6, 4},
    {3, 2, 0, 2, 6, 2, 1, 0, 1, 3, 7, 3, 7, 6, 4, 6, 4, 0, 5, 4, 5, 7, 5, 1},
    {5, 4, 0, 4, 6, 4, 1, 0, 1, 5, 7, 5, 7, 6, 2, 6, 2, 0, 3, 2, 3, 7, 3, 1},
    {6, 2, 0, 2, 3, 2, 4, 0, 4, 6, 7, 6, 7, 3, 1, 3, 1, 0, 5, 1, 5, 7, 5, 4},
    {6, 4, 0, 4, 5, 4, 2, 0, 2, 6, 7, 6, 7, 5, 1, 5, 1, 0, 3, 1, 3, 7, 3, 2},
    {1, 3, 2, 3, 7, 3, 0, 2, 0, 1, 5, 1, 5, 7, 6, 7, 6, 2, 4, 6, 4, 5, 4, 0},
    {1, 5, 4, 5, 7, 5, 0, 4, 0, 1, 3, 1, 3, 7, 6, 7, 6, 4, 2, 6, 2, 3, 2, 0},
    {2, 3, 1, 3, 7, 3, 0, 1, 0, 2, 6, 2, 6, 7, 5, 7, 5, 1, 4, 5, 4, 6, 4, 0},
    {4, 5, 1, 5, 7, 5, 0, 1, 0, 4, 6, 4, 6, 7, 3, 7, 3, 1, 2, 3, 2, 6, 2, 0},
    {2, 6, 4, 6, 7, 6, 0, 4, 0, 2, 3, 2, 3, 7, 5, 7, 5, 4, 1, 5, 1, 3, 1, 0},
    {4, 6, 2, 6, 7, 6, 0, 2, 0, 4, 5, 4, 5, 7, 3, 7, 3, 2, 1, 3, 1, 5, 1, 0},
    {3, 1, 5, 1, 0, 1, 7, 5, 7, 3, 2, 3, 2, 0, 4, 0, 4, 5, 6, 4, 6, 2, 6, 7},
    {5, 1, 3, 1, 0, 1, 7, 3, 7, 5, 4, 5, 4, 0, 2, 0, 2, 3, 6, 2, 6, 4, 6, 7},
    {3, 2, 6, 2, 0, 2, 7, 6, 7, 3, 1, 3, 1, 0, 4, 0, 4, 6, 5, 4, 5, 1, 5, 7},
    {5, 4, 6, 4, 0, 4, 7, 6, 7, 5, 1, 5, 1, 0, 2, 0, 2, 6, 3, 2, 3, 1, 3, 7},
    {6, 2, 3, 2, 0, 2, 7, 3, 7, 6, 4, 6, 4, 0, 1, 0, 1, 3, 5, 1, 5, 4, 5, 7},
    {6, 4, 5, 4, 0, 4, 7, 5, 7, 6, 2, 6, 2, 0, 1, 0, 1, 5, 3, 1, 3, 2, 3, 7},
    {1, 3, 7, 3, 2, 3, 5, 7, 5, 1, 0, 1, 0, 2, 6, 2, 6, 7, 4, 6, 4, 0, 4, 5},
    {1, 5, 7, 5, 4, 5, 3, 7, 3, 1, 0, 1, 0, 4, 6, 4, 6, 7, 2, 6, 2, 0, 2, 3},
    {2, 3, 7, 3, 1, 3, 6, 7, 6, 2, 0, 2, 0, 1, 5, 1, 5, 7, 4, 5, 4, 0, 4, 6},
    {4, 5, 7, 5, 1, 5, 6, 7, 6, 4, 0, 4, 0, 1, 3, 1, 3, 7, 2, 3, 2, 0, 2, 6},
    {2, 6, 7, 6, 4, 6, 3, 7, 3, 2, 0, 2, 0, 4, 5, 4, 5, 7, 1, 5, 1, 0, 1, 3},
    {4, 6, 7, 6, 2, 6, 5, 7, 5, 4, 0, 4, 0, 2, 3, 2, 3, 7, 1, 3, 1, 0, 1, 5},
    {7, 3, 1, 3, 2, 3, 5, 1, 5, 7, 6, 7, 6, 2, 0, 2, 0, 1, 4, 0, 4, 6, 4, 5},
    {7, 5, 1, 5, 4, 5, 3, 1, 3, 7, 6, 7, 6, 4, 0, 4, 0, 1, 2, 0, 2, 6, 2, 3},
    {7, 3, 2, 3, 1, 3, 6, 2, 6, 7, 5, 7, 5, 1, 0, 1, 0, 2, 4, 0, 4, 5, 4, 6},
    {7, 5, 4, 5, 1, 5, 6, 4, 6, 7, 3, 7, 3, 1, 0, 1, 0, 4, 2, 0, 2, 3, 2, 6},
    {7, 6, 2, 6, 4, 6, 3, 2, 3, 7, 5, 7, 5, 4, 0, 4, 0, 2, 1, 0, 1, 5, 1, 3},
    {7, 6, 4, 6, 2, 6, 5, 4, 5, 7, 3, 7, 3, 2, 0, 2, 0, 4, 1, 0, 1, 3, 1, 5},
    {3, 7, 5, 7, 6, 7, 1, 5, 1, 3, 2, 3, 2, 6, 4, 6, 4, 5, 0, 4, 0, 2, 0, 1},
    {5, 7, 3, 7, 6, 7, 1, 3, 1, 5, 4, 5, 4, 6, 2, 6, 2, 3, 0, 2, 0, 4, 0, 1},
    {3, 7, 6, 7, 5, 7, 2, 6, 2, 3, 1, 3, 1, 5, 4, 5, 4, 6, 0, 4, 0, 1, 0, 2},
    {5, 7, 6, 7, 3, 7, 4, 6, 4, 5, 1, 5, 1, 3, 2, 3, 2, 6, 0, 2, 0, 1, 0, 4},
    {6, 7, 3, 7, 5, 7, 2, 3, 2, 6, 4, 6, 4, 5, 1, 5, 1, 3, 0, 1, 0, 4, 0, 2},
    {6, 7, 5, 7, 3, 7, 4, 5, 4, 6, 2, 6, 2, 3, 1, 3, 1, 5, 0, 1, 0, 2, 0, 4},
};

/// @brief Получить индекс прилегающей дочерней ячейки для subface
template <int dim>
constexpr short child(Side<dim> side) {
    if constexpr (dim == 2) { return child_by_subface_2D[side]; }
    else {                    return child_by_subface_3D[side]; }
}

/// @brief Локальные индексы дочерних ячеек, прилегающих к стороне.
template <int dim>
constexpr std::array<short, FpF(dim)> children(Side<dim> side) {
    if constexpr (dim == 2) { return children_by_face_2D[side]; }
    else {                    return children_by_face_3D[side]; }
}

/// @brief Для простой грани и прилегающей дочерней ячейки получить subface
template <int dim>
constexpr Side<dim> subface_by_child(Side<dim> side, int child_idx) {
    if constexpr (dim == 2) { return subface_by_side_child_2D[side][child_idx]; }
    else                    { return subface_by_side_child_3D[side][child_idx]; }
}

/// @brief По subface и повороту получить индекс дочерней соседа
template <int dim>
constexpr int neib_child(Side<dim> side, int symm) {
    if constexpr (dim == 2) { return neib_child_by_symm_subface_2D[symm][side]; }
    else                    { return neib_child_by_symm_subface_3D[symm][side]; }
}

/// @brief По subface и повороту получить индекс дочерней соседа
/// assuming zero rotation
template <int dim>
constexpr int neib_child(Side<dim> side) {
    if constexpr (dim == 2) { return neib_child_by_symm_subface_2D[0][side]; }
    else                    { return neib_child_by_symm_subface_3D[0][side]; }
}

/// @brief По повороту и индексу дочерней (child) определить индекс
/// смежной соседней ячейки через простую сторону
template <int dim>
constexpr int adjacent_child(Side<dim> side, int symm, int child) {
    return neib_child(subface_by_child(side, child), symm);
}

/// @brief По повороту и индексу дочерней (child) определить индекс
/// смежной соседней ячейки через простую сторону
/// assuming zero rotation
template <int dim>
constexpr int adjacent_child(Side<dim> side, int child) {
    return neib_child(subface_by_child(side, child));
}

} // namespace zephyr::geom::indexing