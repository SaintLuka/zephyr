#pragma once

#include <cmath>
#include <iostream>

#include <Eigen/Dense>

namespace zephyr::geom {

typedef Eigen::Matrix<double, 3, 1> Vector3d; ///< 3-ех мерный Eigen вектор-столбец
typedef Eigen::Matrix<double, 4, 1> Vector4d; ///< 4-ех мерный Eigen вектор-столбец
typedef Eigen::Matrix<double, 5, 1> Vector5d; ///< 5-и  мерный Eigen вектор-столбец
typedef Eigen::Matrix<double, 6, 1> Vector6d; ///< 6-и  мерный Eigen вектор-столбец
typedef Eigen::Matrix<double, 7, 1> Vector7d; ///< 7-и  мерный Eigen вектор-столбец

typedef Eigen::Matrix<double, 3, 3> Matrix3d; ///< Диагональная Eigen матрица 3 x 3
typedef Eigen::Matrix<double, 4, 4> Matrix4d; ///< Диагональная Eigen матрица 4 x 4
typedef Eigen::Matrix<double, 5, 5> Matrix5d; ///< Диагональная Eigen матрица 5 x 5
typedef Eigen::Matrix<double, 6, 6> Matrix6d; ///< Диагональная Eigen матрица 6 x 6
typedef Eigen::Matrix<double, 7, 7> Matrix7d; ///< Диагональная Eigen матрица 7 x 7

typedef Eigen::DiagonalMatrix<double, 3, 3> DiagMatrix3d; ///< Диагональная Eigen матрица 3 x 3
typedef Eigen::DiagonalMatrix<double, 4, 4> DiagMatrix4d; ///< Диагональная Eigen матрица 4 x 4
typedef Eigen::DiagonalMatrix<double, 5, 5> DiagMatrix5d; ///< Диагональная Eigen матрица 5 x 5
typedef Eigen::DiagonalMatrix<double, 6, 6> DiagMatrix6d; ///< Диагональная Eigen матрица 6 x 6
typedef Eigen::DiagonalMatrix<double, 7, 7> DiagMatrix7d; ///< Диагональная Eigen матрица 7 x 7

} // namespace zephyr::geom

