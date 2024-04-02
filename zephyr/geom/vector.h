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

/// @brief Установить значения компоненты равными NAN
inline void setNaN(Vector3d& vec) {
    vec.x() = vec.y() = vec.z() = NAN;
}

/// @brief Аналог linspace из numpy.
/// Создает равномерный массив из N чисел на отрезке [t1, t2],
/// концы включаются
inline std::vector<double> linspace(double t1, double t2, size_t N) {
    std::vector<double> res(N);
    double h = (t2 - t1) / (N - 1.0);
    for (size_t i = 0; i < N; ++i) {
        res[i] = t1 + i * h;
    }
    return res;
}

/// @brief Аналог linspace из numpy.
/// Создает равномерный массив из N точек на линии [v1, v2],
/// концы включаются
inline std::vector<Vector3d> linspace(
        const Vector3d& v1, const Vector3d& v2, size_t N) {
    std::vector<Vector3d> res(N);
    Vector3d tau = (v2 - v1) / (N - 1.0);
    for (size_t i = 0; i < N; ++i) {
        res[i] = v1 + tau * i;
    }
    return res;
}

/// @brief Выделить из массива векторов массив x-компонент
inline std::vector<double> get_x(const std::vector<Vector3d>& arr) {
    std::vector<double> res(arr.size());
    for (size_t i = 0; i < arr.size(); ++i) {
        res[i] = arr[i].x();
    }
    return res;
}

/// @brief Выделить из массива векторов массив y-компонент
inline std::vector<double> get_y(const std::vector<Vector3d>& arr) {
    std::vector<double> res(arr.size());
    for (size_t i = 0; i < arr.size(); ++i) {
        res[i] = arr[i].y();
    }
    return res;
}

/// @brief Выделить из массива векторов массив z-компонент
inline std::vector<double> get_z(const std::vector<Vector3d>& arr) {
    std::vector<double> res(arr.size());
    for (size_t i = 0; i < arr.size(); ++i) {
        res[i] = arr[i].z();
    }
    return res;
}

/// @brief Склеить два массива, которые содержат x и y координаты точек,
/// в один массив вершин
inline std::vector<Vector3d> zip(const std::vector<double>& xs,
        const std::vector<double>& ys) {
    if (xs.size() != ys.size()) {
        throw std::runtime_error("zip error: different array size");
    }

    std::vector<Vector3d> res(xs.size());
    for (size_t i = 0; i < xs.size(); ++i) {
        res[i] = {xs[i], ys[i], 0.0};
    }
    return res;
}

/// @brief Склеить три массива, которые содержат (x, y, z) координаты точек,
/// в один массив вершин
inline std::vector<Vector3d> zip(const std::vector<double>& xs,
        const std::vector<double>& ys, const std::vector<double>& zs) {
    if (xs.size() != ys.size() || xs.size() != zs.size()) {
        throw std::runtime_error("zip error: different array size");
    }

    std::vector<Vector3d> res(xs.size());
    for (size_t i = 0; i < xs.size(); ++i) {
        res[i] = {xs[i], ys[i], zs[i]};
    }
    return res;
}

} // namespace zephyr::geom
