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

/// @brief Аналог meshgrid из numpy.
inline std::tuple<std::vector<std::vector<double>>,
                  std::vector<std::vector<double>>> meshgrid(
        const std::vector<double>& x,
        const std::vector<double>& y) {

    std::vector X(x.size(), std::vector<double>(y.size()));
    std::vector Y(x.size(), std::vector<double>(y.size()));

    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = 0; j < y.size(); ++j) {
            X[i][j] = x[i];
            Y[i][j] = y[j];
        }
    }
    return {X, Y};
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
