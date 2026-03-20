// Выставляет повороты на гранях ячеек.
//
// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#include <array>
#include <vector>
#include <format>
#include <iostream>
#include <algorithm>
#include <random>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/geom.h>
#include <zephyr/geom/side.h>
#include <zephyr/geom/indexing.h>
#include <zephyr/mesh/euler/amr_cells.h>

namespace zephyr::mesh::amr {

using namespace geom;

using Vector2i = Eigen::Vector2i;
using Vector3i = Eigen::Vector3i;
using Matrix2i = Eigen::Matrix2i;
using Matrix3i = Eigen::Matrix3i;

// Генератор матриц преобразований для группы симметрии квадрата
inline std::array<Matrix2i, 8> generate_rotation_matrices_2d() {
    std::array<Matrix2i, 8> result{};

    // Вращения (det = +1)
    result[0] = Matrix2i{{ 1,  0}, { 0,  1}}; // 0°
    result[1] = Matrix2i{{ 0, -1}, { 1,  0}}; // 90°
    result[2] = Matrix2i{{-1,  0}, { 0, -1}}; // 180°
    result[3] = Matrix2i{{ 0,  1}, {-1,  0}}; // 270°

    // Отражения (det = -1)
    result[4] = Matrix2i{{ 1,  0}, { 0, -1}}; // отражение по X
    result[5] = Matrix2i{{-1,  0}, { 0,  1}}; // отражение по Y
    result[6] = Matrix2i{{ 0,  1}, { 1,  0}}; // по диагонали
    result[7] = Matrix2i{{ 0, -1}, {-1,  0}}; // по другой диагонали

    return result;
}

// Ключ для сравнения
inline auto matrix_key(const Matrix3i& m) {
    return std::tuple(
        m(0,0), m(0,1), m(0,2),
        m(1,0), m(1,1), m(1,2),
        m(2,0), m(2,1), m(2,2)
    );
};

// Лексикографическое сравнение двух матриц
inline bool compare(const Matrix3i& a, const Matrix3i& b) {
    return matrix_key(a) > matrix_key(b);
}

// Генератор матриц преобразований для группы симметрии куба
inline std::array<Matrix3i, 48> generate_rotation_matrices_3d() {
    std::vector<Matrix3i> rotations;
    std::vector<Matrix3i> reflections;

    std::array perm = {0, 1, 2};
    do {
        for (int sx : {-1, 1})
        for (int sy : {-1, 1})
        for (int sz : {-1, 1}) {
            Matrix3i R = Matrix3i::Zero();

            R(0, perm[0]) = sx;
            R(1, perm[1]) = sy;
            R(2, perm[2]) = sz;

            if (R.determinant() > 0) {
                rotations.push_back(R);
            }
            else {
                reflections.push_back(R);
            }
        }
    } while (std::ranges::next_permutation(perm).found);

    if (rotations.size() != 24 || reflections.size() != 24) {
        throw std::runtime_error("Bad code");
    }

    // Сортировка (тождественная в начале)
    std::ranges::sort(rotations, compare);
    std::ranges::sort(reflections, compare);

    std::array<Matrix3i, 48> result{};

    // сначала вращения
    for (size_t i = 0; i < rotations.size(); ++i) {
        result[i] = rotations[i];
    }

    // затем отражения
    for (size_t i = 0; i < reflections.size(); ++i) {
        result[rotations.size() + i] = reflections[i];
    }

    return result;
}

// Перестановки вершин для группы симметрий квадрата
inline std::array<std::array<int, 4>, 8> generate_permutations_2d() {
    static const std::array vertices_2d = {
        Vector2i{-1, -1},
        Vector2i{+1, -1},
        Vector2i{-1, +1},
        Vector2i{+1, +1}
    };

    // Индексы вершин квадрата [-1, 1]^2 -> 0..4
    auto index = [](Vector2i v) -> int {
        int x = (v.x() + 1) / 2;
        int y = (v.y() + 1) / 2;
        return y * 2 + x;
    };

    auto matrices = generate_rotation_matrices_2d();
    std::array<std::array<int, 4>, 8> perms{};
    for (size_t i = 0; i < matrices.size(); ++i) {
        for (int v = 0; v < 4; ++v) {
            Vector2i transformed = matrices[i] * vertices_2d[v];
            perms[i][v] = index(transformed);
        }
    }
    return perms;
}

// Перестановки вершин для группы симметрий куба
inline std::array<std::array<int, 8>, 48> generate_permutations_3d() {
    static const std::array vertices_3d = {
        Vector3i{-1, -1, -1},
        Vector3i{+1, -1, -1},
        Vector3i{-1, +1, -1},
        Vector3i{+1, +1, -1},
        Vector3i{-1, -1, +1},
        Vector3i{+1, -1, +1},
        Vector3i{-1, +1, +1},
        Vector3i{+1, +1, +1}
    };

    // Индексы вершин куба [-1, 1]^3 -> 0..7
    auto index = [](Vector3i v) -> int {
        int x = (v.x() + 1) / 2;
        int y = (v.y() + 1) / 2;
        int z = (v.z() + 1) / 2;
        return z * 4 + y * 2 + x;
    };

    auto matrices = generate_rotation_matrices_3d();
    std::array<std::array<int, 8>, 48> perms{};
    for (size_t i = 0; i < matrices.size(); ++i) {
        for (int v = 0; v < 8; ++v) {
            Vector3i transformed = matrices[i] * vertices_3d[v];
            perms[i][v] = index(transformed);
        }
    }
    return perms;
}


struct axes2d_t {
    Vector2d ex, ey;

    Matrix2d matrix() const {
        Matrix2d M;
        M.col(0) = ex;
        M.col(1) = ey;
        return M;
    }
};

struct axes3d_t {
    Vector3d ex, ey, ez;

    Matrix3d matrix() const {
        Matrix3d M;
        M.col(0) = ex;
        M.col(1) = ey;
        M.col(2) = ez;
        return M;
    }
};

inline int round_entry(double x) {
    if (x > 0.5) return 1;
    if (x < -0.5) return -1;
    return 0;
}

// Округление к {-1,0,1}
template <int d>
Eigen::Matrix<int, d, d> round_matrix(const Eigen::Matrix<double, d, d>& M) {
    Eigen::Matrix<int, d, d> R;
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
            R(i, j) = round_entry(M(i, j));
        }
    }
    return R;
}

// Поворот квадрата 2 относительно квадрата 1
inline int find_rotation(const axes2d_t& axes1, const axes2d_t& axes2) {
    static const auto group = generate_rotation_matrices_2d();

    Matrix2d A = axes1.matrix();
    Matrix2d B = axes2.matrix();

    // B = R * A => R = B * A^{-1}
    Matrix2d Rf = B * A.inverse();
    Matrix2i R = round_matrix(Rf);

    // ищем в группе
    auto it = std::ranges::find(group, R);
    if (it != group.end()) {
        return std::distance(group.begin(), it);
    }
    return -1;
}

// Выполнить перестановку в контейнере
template <typename T>
T permute(const T& container, std::span<const int> perm) {
    T result{};
    for (int i = 0; i < perm.size(); ++i) {
        result[i] = container[perm[i]];
    }
    return result;
}

// Выполнить обратную перестановку в контейнере
template <typename T>
T permute_inv(const T& container, std::span<const int> perm) {
    T result{};
    for (int i = 0; i < perm.size(); ++i) {
        result[perm[i]] = container[i];
    }
    return result;
}

// Выполнить перестановку в контейнере
template <typename T>
T permute_2d(const T& container, int r) {
    static const auto group = generate_permutations_2d();
    return permute(container, group[r]);
}

// Выполнить обратную перестановку в контейнере
template <typename T>
T permute_inv_2d(const T& container, int r) {
    static const auto group = generate_permutations_2d();
    return permute_inv(container, group[r]);
}

// Выполнить перестановку в контейнере
template <typename T>
T permute_3d(const T& container, int r) {
    static const auto group = generate_permutations_3d();
    return permute(container, group[r]);
}

// Выполнить перестановку в контейнере
template <typename T>
T permute_inv_3d(const T& container, int r) {
    static const auto group = generate_permutations_3d();
    return permute_inv(container, group[r]);
}

// Поворот куба 2 относительно куба 1
inline int find_rotation(const axes3d_t& axes1, const axes3d_t& axes2) {
    static const auto group = generate_rotation_matrices_3d();

    Matrix3d A = axes1.matrix();
    Matrix3d B = axes2.matrix();

    // B = R * A => R = B * A^{-1}
    Matrix3d Rf = B * A.inverse();
    Matrix3i R = round_matrix(Rf);

    // ищем в группе
    auto it = std::ranges::find(group, R);
    if (it != group.end()) {
        return std::distance(group.begin(), it);
    }
    return -1;
}

inline std::uint8_t find_rotation(const Quad& quad1, Side2D side, const Quad& quad2) {
    // Индексы вершин квадрата [-1, 1]^2 -> 0..4
    auto vs = [](int i, int j) -> int {
        int x = (i + 1) / 2;
        int y = (j + 1) / 2;
        return y * 2 + x;
    };
    
    constexpr int face_nodes[4][2] = {
        {vs(-1, -1), vs(-1, +1)},
        {vs(+1, -1), vs(+1, +1)},
        {vs(-1, -1), vs(+1, -1)},
        {vs(-1, +1), vs(+1, +1)},
    };

    static const auto permutations_2d = generate_permutations_2d();

    Side2D twin = (side % 2 == 0) ? side + 1 : side - 1;    

    const Vector3d& v1 = quad1[face_nodes[side][0]];
    const Vector3d& v2 = quad1[face_nodes[side][1]];

    for (int r = 0; r < permutations_2d.size(); ++r) {
        Quad quad3 = permute_inv_2d(quad2, r);

        const Vector3d& w1 = quad3[face_nodes[twin][0]];
        const Vector3d& w2 = quad3[face_nodes[twin][1]];
        if ((v1 - w1).norm() < 1.0e-13 && (v2 - w2).norm() < 1.0e-13) {
            return r;
        }
    }
    return -1;
}

inline std::uint8_t find_rotation(const Cube& cube1, Side3D side, const Cube& cube2) {
    // Индексы вершин куба [-1, 1]^3 -> 0..7
    auto vs = [](int i, int j, int k) -> int {
        int x = (i + 1) / 2;
        int y = (j + 1) / 2;
        int z = (k + 1) / 2;
        return z * 4 + y * 2 + x;
    };
    
    constexpr int face_nodes[6][4] = {
        {vs(-1, -1, -1), vs(-1, +1, -1), vs(-1, -1, +1), vs(-1, +1, +1)},
        {vs(+1, -1, -1), vs(+1, +1, -1), vs(+1, -1, +1), vs(+1, +1, +1)},
        {vs(-1, -1, -1), vs(+1, -1, -1), vs(-1, -1, +1), vs(+1, -1, +1)},
        {vs(-1, +1, -1), vs(+1, +1, -1), vs(-1, +1, +1), vs(+1, +1, +1)},
        {vs(-1, -1, -1), vs(+1, -1, -1), vs(-1, +1, -1), vs(+1, +1, -1)},
        {vs(-1, -1, +1), vs(+1, -1, +1), vs(-1, +1, +1), vs(+1, +1, +1)},
    };

    static const auto permutations_3d = generate_permutations_3d();

    Side2D twin = (side % 2 == 0) ? side + 1 : side - 1;

    const Vector3d& v1 = cube1[face_nodes[side][0]];
    const Vector3d& v2 = cube1[face_nodes[side][1]];
    const Vector3d& v3 = cube1[face_nodes[side][2]];

    for (int r = 0; r < permutations_3d.size(); ++r) {
        Cube cube3 = permute_inv_3d(cube2, r);

        const Vector3d& w1 = cube3[face_nodes[twin][0]];
        const Vector3d& w2 = cube3[face_nodes[twin][1]];
        const Vector3d& w3 = cube3[face_nodes[twin][2]];
        if ((v1 - w1).norm() < 1.0e-13 &&
            (v2 - w2).norm() < 1.0e-13 &&
            (v3 - w3).norm() < 1.0e-13) {
            return r;
        }
    }
    return -1;
}

template <int dim>
using Map = std::conditional_t<dim == 2, Quad, Cube>;

template <int dim>
void find_rotations_impl(AmrCells& locals, const AmrCells& aliens) {
    z_assert(locals.adaptive(), "find_rotations: not adaptive mesh");
    for (index_t ic = 0; ic < locals.size(); ++ic) {
        Map<dim> map1 = locals.mapping<dim>(ic).reduce();

        for (Side<dim> side: Side<dim>::items()) {
            index_t iface = locals.face_begin[ic] + side;
            if (locals.faces.is_boundary(iface)) {
                locals.faces.adjacent.rotation[iface] = 0;
                continue;
            }

            auto [neibs, jc] = locals.faces.adjacent.get_neib(iface, locals, aliens);

            Map<dim> map2 = neibs.mapping<dim>(jc).reduce();

            auto r = find_rotation(map1, side, map2);
            if (r == 255) {
                throw std::runtime_error("Invalid rotation");
            }
            locals.faces.adjacent.rotation[iface] = r;
        }
    }
}

/// @brief Автоматический выбор размерности
inline void find_rotations(AmrCells &locals) {
    static AmrCells aliens;

    if (locals.empty()) return;

    if (locals.dim() < 3) {
        amr::find_rotations_impl<2>(locals, aliens);
    }
    else {
        amr::find_rotations_impl<3>(locals, aliens);
    }
}

} // namespace zephyr::mesh::amr