#pragma once

#include <cmath>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include <zephyr/configuration.h>

#ifndef ZEPHYR_EIGEN
#include <zephyr/geom/eigen/interface.h>
#endif

/// @brief Геометрические примитивы, поиск площадей, объемов, сечений.
namespace zephyr::geom {

#ifdef ZEPHYR_EIGEN
using Vector3d = Eigen::Matrix<double, 3, 1>; ///< 3-ех мерный Eigen вектор-столбец
using Matrix3d = Eigen::Matrix<double, 3, 3>; ///< Eigen матрица 4 x 4
#else
using Vector3d = eigen_wrapper::Matrix<double, 3, 1>;
using Matrix3d = eigen_wrapper::Matrix<double, 3, 3>;
#endif

using Vector2d = Eigen::Matrix<double, 2, 1>; ///< 2-х мерный Eigen вектор-столбец
using Vector4d = Eigen::Matrix<double, 4, 1>; ///< 4-ех мерный Eigen вектор-столбец
using Vector5d = Eigen::Matrix<double, 5, 1>; ///< 5-и  мерный Eigen вектор-столбец
using Vector6d = Eigen::Matrix<double, 6, 1>; ///< 6-и  мерный Eigen вектор-столбец
using Vector7d = Eigen::Matrix<double, 7, 1>; ///< 7-и  мерный Eigen вектор-столбец

using Matrix2d = Eigen::Matrix<double, 4, 4>; ///< Eigen матрица 2 x 2
using Matrix4d = Eigen::Matrix<double, 4, 4>; ///< Eigen матрица 4 x 4
using Matrix5d = Eigen::Matrix<double, 5, 5>; ///< Eigen матрица 5 x 5
using Matrix6d = Eigen::Matrix<double, 6, 6>; ///< Eigen матрица 6 x 6
using Matrix7d = Eigen::Matrix<double, 7, 7>; ///< Eigen матрица 7 x 7

using DiagMatrix3d = Eigen::DiagonalMatrix<double, 3, 3>; ///< Диагональная Eigen матрица 3 x 3
using DiagMatrix4d = Eigen::DiagonalMatrix<double, 4, 4>; ///< Диагональная Eigen матрица 4 x 4
using DiagMatrix5d = Eigen::DiagonalMatrix<double, 5, 5>; ///< Диагональная Eigen матрица 5 x 5
using DiagMatrix6d = Eigen::DiagonalMatrix<double, 6, 6>; ///< Диагональная Eigen матрица 6 x 6
using DiagMatrix7d = Eigen::DiagonalMatrix<double, 7, 7>; ///< Диагональная Eigen матрица 7 x 7

/// @brief Eigen вектор-столбец
template <typename Scalar, int Rows>
using Vector = Eigen::Matrix<Scalar, Rows, 1>;

/// @brief Eigen вектор-столбец
template <typename Scalar, int Rows>
using Array = Eigen::Array<Scalar, Rows, 1>;

/// @brief Eigen матрица
template <typename Scalar, int Rows, int Cols>
using Matrix = Eigen::Matrix<Scalar, Rows, Cols>;


/// @brief Тип произвольного класса, приведенный к Eigen::Array
template<class T>
using ei_arr = Eigen::Array<double, sizeof(T) / sizeof(double), 1>;

/// @brief Тип произвольного класса, приведенный к Eigen::Matrix
template<class T>
using ei_vec = Eigen::Matrix<double, sizeof(T) / sizeof(double), 1>;



} // namespace zephyr::geom
