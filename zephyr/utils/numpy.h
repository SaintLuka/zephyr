#pragma once

#include <numeric>
#include <vector>
#include <zephyr/geom/vector.h>

/// @brief Несколько удобных функций numpy, которые работают с `std::vector`.
namespace zephyr::np {

// Типы для обращения np::real, np::vec3
using real = double;
using vec3 = geom::Vector3d;

template <typename T>
using Array = std::vector<T>;

using Array1d = std::vector<real>;
using Array1v = std::vector<vec3>;

using Array2d = std::vector<std::vector<real>>;


/// @{ @name Инициализация массивов

/// @brief Массив длины N, заполненный нулями
template <typename T = np::real>
Array<T> zeros(size_t N);

template <>
Array<real> zeros<real>(size_t N);

template <>
Array<vec3> zeros<vec3>(size_t N);


/// @brief Двумерный массив из нулей длины N x M
Array2d zeros(size_t N, size_t M);

Array1d zeros_like(const Array1d& arr);

Array2d zeros_like(const Array2d& arr);

template <typename T>
Array<T> empty_like(const Array<T>& arr) { return Array<T>(arr.size()); }

inline Array2d empty_like(const Array2d& arr) {
    if (arr.empty()) { return {}; }
    return Array2d(arr.size(), Array1d(arr[0].size()));
}

/// @}

/// @{ @name Операции свёртки

template <typename T>
T min(const std::vector<T>& vec) { return *std::min_element(vec.begin(), vec.end()); }

template <typename T>
T max(const std::vector<T>& vec) { return *std::max_element(vec.begin(), vec.end()); }

template <typename T>
std::enable_if_t<std::is_arithmetic_v<T>, T>
mean(const std::vector<T>& vec) {
    T sum = std::accumulate(vec.begin(), vec.end(), T{});
    return sum / vec.size();
}

/// @}

/// @{ @name Генерация сеток

/// @brief
Array<real> linspace(real t1, real t2, size_t N);

/// @brief Аналог linspace из numpy.
/// Создает равномерный массив из N точек на линии [v1, v2],
/// концы включаются
Array<vec3> linspace(const vec3& v1, const vec3& v2, size_t N);

/// @brief
std::tuple<Array2d, Array2d> meshgrid(const Array1d& x, const Array1d& y);

/// @}

/// @{ @name Объединение массивов и выделение осей

/// @brief Склеить два массива, которые содержат x и y координаты точек,
/// в один массив вершин
Array1v zip(const Array1d& xs, const Array1d& ys);

/// @brief Склеить три массива, которые содержат (x, y, z) координаты точек,
/// в один массив вершин
Array1v zip(const Array1d& xs, const Array1d& ys, const Array1d& zs);

/// @brief Выделить из массива векторов массив x-компонент
Array1d get_x(const Array1v& arr);

/// @brief Выделить из массива векторов массив y-компонент
Array1d get_y(const Array1v& arr);

/// @brief Выделить из массива векторов массив z-компонент
Array1d get_z(const Array1v& arr);

/// @}

} // namespace zephyr::np


