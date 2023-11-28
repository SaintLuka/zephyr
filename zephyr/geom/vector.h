#pragma once

#include <cmath>

#include <zephyr/configuration.h>

#ifdef ZEPHYR_ENABLE_EIGEN
#include <Eigen/Dense>
#else
#include <zephyr/geom/eigen/interface.h>
#endif

namespace zephyr::geom {

#ifdef ZEPHYR_ENABLE_EIGEN
using namespace Eigen;
#endif

typedef Matrix<double, 3, 1> Vector3d; ///< 3-ех мерный Eigen вектор-столбец
typedef Matrix<double, 4, 1> Vector4d; ///< 4-ех мерный Eigen вектор-столбец
typedef Matrix<double, 5, 1> Vector5d; ///< 5-и  мерный Eigen вектор-столбец
typedef Matrix<double, 6, 1> Vector6d; ///< 6-и  мерный Eigen вектор-столбец
typedef Matrix<double, 7, 1> Vector7d; ///< 7-и  мерный Eigen вектор-столбец

typedef Matrix<double, 3, 3> Matrix3d; ///< Eigen матрица 3 x 3
typedef Matrix<double, 4, 4> Matrix4d; ///< Eigen матрица 4 x 4
typedef Matrix<double, 5, 5> Matrix5d; ///< Eigen матрица 5 x 5
typedef Matrix<double, 6, 6> Matrix6d; ///< Eigen матрица 6 x 6
typedef Matrix<double, 7, 7> Matrix7d; ///< Eigen матрица 7 x 7

typedef DiagonalMatrix<double, 3, 3> DiagMatrix3d; ///< Диагональная Eigen матрица 3 x 3
typedef DiagonalMatrix<double, 4, 4> DiagMatrix4d; ///< Диагональная Eigen матрица 4 x 4
typedef DiagonalMatrix<double, 5, 5> DiagMatrix5d; ///< Диагональная Eigen матрица 5 x 5
typedef DiagonalMatrix<double, 6, 6> DiagMatrix6d; ///< Диагональная Eigen матрица 6 x 6
typedef DiagonalMatrix<double, 7, 7> DiagMatrix7d; ///< Диагональная Eigen матрица 7 x 7

inline void setNaN(Vector3d& vec) {
    vec.x() = vec.y() = vec.z() = NAN;
}

} // namespace zephyr::geom
