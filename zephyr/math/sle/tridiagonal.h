/// @file Решение систем линейных уравнений с трехдиагональными матрицами
/// x = solve(A, B, C, F) -- аргументы не меняются, решение возвращается.
/// solve_l(A, B, C, F)   -- решение записывается в последний массив F.
/// solve_m(A, B, C, F)   -- аргументы изменяются, решение записывается
/// в последний массив F.
///
/// solve_cyclic -- решения для почти трехдиагональных матриц со связью
/// первого и последнего элемента, возникают при периодических граничных
/// условиях. Наименования аналогичные.
/// x = solve_cyclic(A, B, C, F) -- аргументы не меняются, решение возвращается.
/// solve_cyclic_l(A, B, C, F)   -- решение записывается в последний массив F.
/// solve_cyclic_m(A, B, C, F)   -- аргументы изменяются, решение записывается
/// в последний массив F.

#pragma once

#include <vector>

namespace zephyr::math::tridiagonal {

using array = std::vector<double>;

/// @brief Решить систему с трехдиагональной матрицей.
/// B_0 x_0 + C_0 x_1 = F_0
/// A_i x_{i - 1} + B_i x_i + C_i x_{i + 1} = F_i, i = 1..n-1
/// A_n x_{n - 1} + B_n x_n = F_n
array solve(const array& A, const array& B, const array& C, const array& F);

/// @brief Решить систему с трехдиагональной матрицей.
/// Решение записывается в последний массив F.
/// B_0 x_0 + C_0 x_1 = F_0
/// A_i x_{i - 1} + B_i x_i + C_i x_{i + 1} = F_i, i = 1..n-1
/// A_n x_{n - 1} + B_n x_n = F_n
void solve_l(const array& A, const array& B, const array& C, array& F);

/// @brief Решить систему с трехдиагональной матрицей, изменяются последние
/// два аргумента (массивы C и F), решение записывается в последний массив F.
/// B_0 x_0 + C_0 x_1 = F_0
/// A_i x_{i - 1} + B_i x_i + C_i x_{i + 1} = F_i, i = 1..n-1
/// A_n x_{n - 1} + B_n x_n = F_n
void solve_m(const array& A, const array& B, array& C, array& F);

/// @brief Решить систему с почти трехдиагональной матрицей и периодическими
/// граничными условиями.
/// A_0 x_n + B_0 x_0 + C_0 x_1 = F_0
/// A_i x_{i - 1} + B_i x_i + C_i x_{i + 1} = F_i, i = 1..n-1
/// A_n x_{n - 1} + B_n x_n + C_n x_0 = F_n
array solve_cyclic(const array& A, const array& B, const array& C, const array& F);

/// @brief Решить систему с почти трехдиагональной матрицей и периодическими
/// граничными условиями. Решение записывается в последний массив F.
/// A_0 x_n + B_0 x_0 + C_0 x_1 = F_0
/// A_i x_{i - 1} + B_i x_i + C_i x_{i + 1} = F_i, i = 1..n-1
/// A_n x_{n - 1} + B_n x_n + C_n x_0 = F_n
void solve_cyclic_l(const array& A, const array& B, const array& C, array& F);

/// @brief Решить систему с почти трехдиагональной матрицей и периодическими
/// граничными условиями, изменяются последние три аргумента (массивы B, C и F),
/// решение записывается в последний массив F.
/// A_0 x_n + B_0 x_0 + C_0 x_1 = F_0
/// A_i x_{i - 1} + B_i x_i + C_i x_{i + 1} = F_i, i = 1..n-1
/// A_n x_{n - 1} + B_n x_n + C_n x_0 = F_n
void solve_cyclic_m(const array& A, array& B, array& C, array& F);

} // namespace zephyr::math::tridiagonal