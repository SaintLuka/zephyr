#pragma once

#include <array>
#include <string>

namespace zephyr::mesh {

/// @brief Число вершин простой грани.
constexpr int VpF(int dim) { return dim < 3 ? 2 : 4; }

/// @brief Количество подграней.
constexpr int FpF(int dim) { return dim < 3 ? 2 : 4; }

/// @brief Число граней у простой ячейки.
constexpr int FpC(int dim) { return dim < 3 ? 4 : 6; }

/// @brief Количество дочерних ячеек.
constexpr int CpC(int dim) { return dim < 3 ? 4 : 8; }

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

/// @brief Индексация вершин в кубической AMR-ячейке. Повторяет индексацию функции SqCube::iss
template<short i, short j, short k = -1>
static constexpr short iss() {
    static_assert(i * i <= 1 && j * j <= 1 && k * k <= 1, "Available indices: {-1, 0, 1}");
    return 9 * (k + 1) + 3 * (j + 1) + i + 1;
}

static constexpr std::array<short, 8> simple_face_indices_2D[4] = {
    {iss<-1, -1>(), iss<-1, +1>(), -1, -1, -1, -1, -1, -1},
    {iss<+1, -1>(), iss<+1, +1>(), -1, -1, -1, -1, -1, -1},
    {iss<-1, -1>(), iss<+1, -1>(), -1, -1, -1, -1, -1, -1},
    {iss<-1, +1>(), iss<+1, +1>(), -1, -1, -1, -1, -1, -1},
};

static constexpr std::array<short, 8> simple_face_indices_3D[6] = {
    {iss<-1, -1, -1>(), iss<-1, +1, -1>(), iss<-1, -1, +1>(), iss<-1, +1, +1>(), -1, -1, -1, -1},
    {iss<+1, -1, -1>(), iss<+1, +1, -1>(), iss<+1, -1, +1>(), iss<+1, +1, +1>(), -1, -1, -1, -1},
    {iss<-1, -1, -1>(), iss<+1, -1, -1>(), iss<-1, -1, +1>(), iss<+1, -1, +1>(), -1, -1, -1, -1},
    {iss<-1, +1, -1>(), iss<+1, +1, -1>(), iss<-1, +1, +1>(), iss<+1, +1, +1>(), -1, -1, -1, -1},
    {iss<-1, -1, -1>(), iss<+1, -1, -1>(), iss<-1, +1, -1>(), iss<+1, +1, -1>(), -1, -1, -1, -1},
    {iss<-1, -1, +1>(), iss<+1, -1, +1>(), iss<-1, +1, +1>(), iss<+1, +1, +1>(), -1, -1, -1, -1},
};

static constexpr std::array<short, 8> complex_face_indices_2D[8] = {
    {iss<-1,-1>(), iss<-1, 0>(), -1, -1, -1, -1, -1, -1}, // L[0]
    {iss<+1,-1>(), iss<+1, 0>(), -1, -1, -1, -1, -1, -1}, // R[0]
    {iss<-1,-1>(), iss< 0,-1>(), -1, -1, -1, -1, -1, -1}, // B[0]
    {iss<-1,+1>(), iss< 0,+1>(), -1, -1, -1, -1, -1, -1}, // T[0]
                    
    {iss<-1, 0>(), iss<-1,+1>(), -1, -1, -1, -1, -1, -1}, // L[1]
    {iss<+1, 0>(), iss<+1,+1>(), -1, -1, -1, -1, -1, -1}, // R[1]
    {iss< 0,-1>(), iss<+1,-1>(), -1, -1, -1, -1, -1, -1}, // B[1]
    {iss< 0,+1>(), iss<+1,+1>(), -1, -1, -1, -1, -1, -1}, // T[1]
};

static constexpr std::array<short, 8> complex_face_indices_3D[24] = {
    {iss<-1,-1,-1>(), iss<-1, 0,-1>(), iss<-1,-1, 0>(), iss<-1, 0, 0>(), -1, -1, -1, -1}, // L[0]
    {iss<+1,-1,-1>(), iss<+1, 0,-1>(), iss<+1,-1, 0>(), iss<+1, 0, 0>(), -1, -1, -1, -1}, // R[0]
    {iss<-1,-1,-1>(), iss< 0,-1,-1>(), iss<-1,-1, 0>(), iss< 0,-1, 0>(), -1, -1, -1, -1}, // B[0]
    {iss<-1,+1,-1>(), iss< 0,+1,-1>(), iss<-1,+1, 0>(), iss< 0,+1, 0>(), -1, -1, -1, -1}, // T[0]
    {iss<-1,-1,-1>(), iss<0, -1,-1>(), iss<-1, 0,-1>(), iss< 0, 0,-1>(), -1, -1, -1, -1}, // Z[0]
    {iss<-1,-1,+1>(), iss< 0,-1,+1>(), iss<-1, 0,+1>(), iss< 0, 0,+1>(), -1, -1, -1, -1}, // F[0]
                    
    {iss<-1, 0,-1>(), iss<-1,+1,-1>(), iss<-1, 0, 0>(), iss<-1,+1, 0>(), -1, -1, -1, -1}, // L[1]
    {iss<+1, 0,-1>(), iss<+1,+1,-1>(), iss<+1, 0, 0>(), iss<+1,+1, 0>(), -1, -1, -1, -1}, // R[1]
    {iss< 0,-1,-1>(), iss<+1,-1,-1>(), iss< 0,-1, 0>(), iss<+1,-1, 0>(), -1, -1, -1, -1}, // B[1]
    {iss< 0,+1,-1>(), iss<+1,+1,-1>(), iss< 0,+1, 0>(), iss<+1,+1, 0>(), -1, -1, -1, -1}, // T[1]
    {iss< 0,-1,-1>(), iss<+1,-1,-1>(), iss< 0, 0,-1>(), iss<+1, 0,-1>(), -1, -1, -1, -1}, // Z[1]
    {iss< 0,-1,+1>(), iss<+1,-1,+1>(), iss< 0, 0,+1>(), iss<+1, 0,+1>(), -1, -1, -1, -1}, // F[1]
                    
    {iss<-1,-1, 0>(), iss<-1, 0, 0>(), iss<-1,-1,+1>(), iss<-1, 0,+1>(), -1, -1, -1, -1}, // L[2]
    {iss<+1,-1, 0>(), iss<+1, 0, 0>(), iss<+1,-1,+1>(), iss<+1, 0,+1>(), -1, -1, -1, -1}, // R[2]
    {iss<-1,-1, 0>(), iss< 0,-1, 0>(), iss<-1,-1,+1>(), iss< 0,-1,+1>(), -1, -1, -1, -1}, // B[2]
    {iss<-1,+1, 0>(), iss< 0,+1, 0>(), iss<-1,+1,+1>(), iss< 0,+1,+1>(), -1, -1, -1, -1}, // T[2]
    {iss<-1, 0,-1>(), iss< 0, 0,-1>(), iss<-1,+1,-1>(), iss< 0,+1,-1>(), -1, -1, -1, -1}, // Z[2]
    {iss<-1, 0,+1>(), iss< 0, 0,+1>(), iss<-1,+1,+1>(), iss< 0,+1,+1>(), -1, -1, -1, -1}, // F[2]
                    
    {iss<-1, 0, 0>(), iss<-1,+1, 0>(), iss<-1, 0,+1>(), iss<-1,+1,+1>(), -1, -1, -1, -1}, // L[3]
    {iss<+1, 0, 0>(), iss<+1,+1, 0>(), iss<+1, 0,+1>(), iss<+1,+1,+1>(), -1, -1, -1, -1}, // R[3]
    {iss< 0,-1, 0>(), iss<+1,-1, 0>(), iss< 0,-1,+1>(), iss<+1,-1,+1>(), -1, -1, -1, -1}, // B[3]
    {iss< 0,+1, 0>(), iss<+1,+1, 0>(), iss< 0,+1,+1>(), iss<+1,+1,+1>(), -1, -1, -1, -1}, // T[3]
    {iss< 0, 0,-1>(), iss<+1, 0,-1>(), iss< 0,+1,-1>(), iss<+1,+1,-1>(), -1, -1, -1, -1}, // Z[3]
    {iss< 0, 0,+1>(), iss<+1, 0,+1>(), iss< 0,+1,+1>(), iss<+1,+1,+1>(), -1, -1, -1, -1}, // F[3]
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

/// @brief Структура позволяет хранить индексы граней AMR-ячейки.
/// В большинстве контекстов работает как простой `enum Side : int`,
/// то есть позволяет получить `Side2D::LEFT`, `Side3D::BACK`, ...
/// Также позволяет получить индексы подграней, к примеру,
/// для трехмерной ячейки `Side3D::BACK[0]`, `Side3D::BACK[1]`, ...
/// @tparam dim Размерность ячейки
template <int dim>
struct Side final {
    // Только двумерные и трехмерные
    static_assert(dim == 2 || dim == 3);

    /// @brief Именованные индексы простых граней
    static const Side LEFT;    //=0
    static const Side RIGHT;   //=1
    static const Side BOTTOM;  //=2
    static const Side TOP;     //=3
    static const Side BACK;    //=4; Только для dim == 3.
    static const Side FRONT;   //=5; Только для dim == 3.

    /// @brief Сокращенные имена граней (Z = BACK)
    static const Side L, R, B, T, Z, F;

    /// @brief Индекс грани
    int val;

    /// @brief Неявное создание из типа int
    constexpr Side(int v) noexcept : val(v) { }

    /// @brief Неявное преобразование в int
    constexpr operator int() const { return val; }

    /// @brief Число простых граней ячейки
    static constexpr int count() { return 2 * dim; }

    /// @brief Полное число подграней ячейки
    static constexpr int n_subfaces() { return count() * FpF(dim); }

    /// @brief Индексы подграней на стороне в списке faces.
    std::array<Side, FpF(dim)> subfaces() const {
        if constexpr (dim == 2) {
            return {Side{val}, Side{val}[1]};
        } else {
            return {Side{val}, Side{val}[1], Side{val}[2], Side{val}[3]};
        }
    }

    /// @brief Получить подгрань сложной грани
    constexpr Side operator[](int idx) const {
        return Side(val + count() * idx);
    }

    /// @brief Перечисление простых граней двумерной ячейки
    static constexpr std::array<Side, count()> items() {
        if constexpr (dim == 2) {
            return {LEFT, RIGHT, BOTTOM, TOP};
        } else {
            return {LEFT, RIGHT, BOTTOM, TOP, BACK, FRONT};
        }
    }

    /// @brief Простая сторона из возможной подграни
    constexpr Side simple() const { return val % count(); }

    /// @brief Получить сторону по нормальному вектору направления
    /// Side3D::LEFT == Side3D::by_dir<-1, 0, 0>().
    template <int i, int j, int k = 0>
    static constexpr Side by_dir() {
        if constexpr (i < 0 && j == 0 && k == 0) { return LEFT; }
        if constexpr (i > 0 && j == 0 && k == 0) { return RIGHT; }
        if constexpr (i == 0 && j < 0 && k == 0) { return BOTTOM; }
        if constexpr (i == 0 && j > 0 && k == 0) { return TOP; }
        if constexpr (i == 0 && j == 0 && k < 0) { return BACK; }
        if constexpr (i == 0 && j == 0 && k > 0) { return FRONT; }
        return -1;
    }

    /// @brief Преобразование в строку
    std::string to_string() const {
        std::string num = "[" + std::to_string(val / count()) + "]";
        switch (simple()) {
            case 0: return "LEFT" + num;
            case 1: return "RIGHT" + num;
            case 2: return "BOTTOM" + num;
            case 3: return "TOP" + num;
            case 4: return "BACK" + num;
            case 5: return "FRONT" + num;
            default: return "STRANGE FACE";
        }
    }

    // Для AMR-ячейки вершины и грани упорядочены определенным образом,
    // поэтому индексы вершин на гранях заведомо известны.

    /// @brief Индексы вершин простой грани в AMR-ячейке
    /// sf -- Simple Face
    constexpr std::array<short, 8> sf() const {
        if constexpr (dim == 2) { return simple_face_indices_2D[val]; }
        else                    { return simple_face_indices_3D[val]; }
    }

    /// @brief Индексы вершин части составной грани AMR-ячейки
    /// cf -- Complex Face
    constexpr std::array<short, 8> cf() const {
        if constexpr (dim == 2) { return complex_face_indices_2D[val]; }
        else                    { return complex_face_indices_3D[val]; }
    }

    /// @brief Получить индекс прилегающей дочерней ячейки для subface
    constexpr short child() const {
        if constexpr (dim == 2) { return child_by_subface_2D[val]; }
        else {                    return child_by_subface_3D[val]; }
    }

    /// @brief Локальные индексы дочерних ячеек, прилегающих к стороне.
    constexpr std::array<short, FpF(dim)> children() {
        if constexpr (dim == 2) { return children_by_face_2D[val]; }
        else {                    return children_by_face_3D[val]; }
    }

    /// @brief Для простой грани и прилегающей дочерней ячейки получить subface
    constexpr Side subface_by_child(int child_idx) const {
        if constexpr (dim == 2) { return subface_by_side_child_2D[val][child_idx]; }
        else                    { return subface_by_side_child_3D[val][child_idx]; }
    }

    /// @brief По subface и повороту получить индекс дочерней соседа
    constexpr int neib_child(int symm) const {
        if constexpr (dim == 2) { return neib_child_by_symm_subface_2D[symm][val]; }
        else                    { return neib_child_by_symm_subface_3D[symm][val]; }
    }

    /// @brief По subface и повороту получить индекс дочерней соседа
    /// assuming zero rotation
    constexpr int neib_child() const {
        if constexpr (dim == 2) { return neib_child_by_symm_subface_2D[0][val]; }
        else                    { return neib_child_by_symm_subface_3D[0][val]; }
    }

    /// @brief По повороту и индексу дочерней (child) определить индекс
    /// смежной соседней ячейки через простую сторону
    constexpr int adjacent_child(int symm, int child) const {
        return subface_by_child(child).neib_child(symm);
    }

    /// @brief По повороту и индексу дочерней (child) определить индекс
    /// смежной соседней ячейки через простую сторону
    /// assuming zero rotation
    constexpr int adjacent_child(int child) const {
        return subface_by_child(child).neib_child();
    }
};

// Инициализация констант для Side<2>, ошибка компиляции
// при попытке использовать Side2D::BACK
using Side2D = Side<2>;
template <> constexpr Side2D Side2D::LEFT   = 0;
template <> constexpr Side2D Side2D::RIGHT  = 1;
template <> constexpr Side2D Side2D::BOTTOM = 2;
template <> constexpr Side2D Side2D::TOP    = 3;
template <> constexpr Side2D Side2D::L = LEFT;
template <> constexpr Side2D Side2D::R = RIGHT;
template <> constexpr Side2D Side2D::B = BOTTOM;
template <> constexpr Side2D Side2D::T = TOP;

// Инициализация констант для Side<3>
using Side3D = Side<3>;
template <> constexpr Side3D Side3D::LEFT   = 0;
template <> constexpr Side3D Side3D::RIGHT  = 1;
template <> constexpr Side3D Side3D::BOTTOM = 2;
template <> constexpr Side3D Side3D::TOP    = 3;
template <> constexpr Side3D Side3D::BACK   = 4;
template <> constexpr Side3D Side3D::FRONT  = 5;
template <> constexpr Side3D Side3D::L = LEFT;
template <> constexpr Side3D Side3D::R = RIGHT;
template <> constexpr Side3D Side3D::B = BOTTOM;
template <> constexpr Side3D Side3D::T = TOP;
template <> constexpr Side3D Side3D::Z = BACK;
template <> constexpr Side3D Side3D::F = FRONT;

/// @brief Вывод названия грани в поток
template <int dim>
std::ostream& operator<<(std::ostream& os, Side<dim> side) {
    os << side.to_string();
    return os;
}

/// @brief Индекс грани в строку, dim -- размерность
inline std::string side_to_string(int s, int dim) {
    return dim == 2 ? Side2D(s).to_string() : Side3D(s).to_string();
}

/// @brief Требуемый порядок граней
static_assert(Side2D::LEFT   == 0, "Left   != 0");
static_assert(Side2D::RIGHT  == 1, "Right  != 1");
static_assert(Side2D::BOTTOM == 2, "Bottom != 2");
static_assert(Side2D::TOP    == 3, "Top    != 3");

static_assert(Side3D::LEFT   == 0, "Left   != 0");
static_assert(Side3D::RIGHT  == 1, "Right  != 1");
static_assert(Side3D::BOTTOM == 2, "Bottom != 2");
static_assert(Side3D::TOP    == 3, "Top    != 3");
static_assert(Side3D::BACK   == 4, "Back   != 4");
static_assert(Side3D::FRONT  == 5, "Front  != 5");

} // namespace zephyr::mesh