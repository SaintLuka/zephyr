#pragma once

#include <array>
#include <vector>

#include <zephyr/geom/vector.h>

namespace zephyr::geom {

/// @brief Список вершин произвольной длины
/// для описания многоугольников
using VerticesList = std::vector<Vector3d>;

/// @brief Площадь произвольного многоугольника
/// Вершины должны располагаться в плоскости (x, y).
double area(const VerticesList& vs);

/// @brief Центр произвольного многоугольника (среднее вершин)
Vector3d center(const VerticesList& vs);

/// @brief Барицентр произвольного многоугольника.
/// Вершины должны располагаться в плоскости (x, y).
/// @param area Площадь многоугольника
Vector3d centroid(const VerticesList &vs, double area = 0.0);

} // namespace zephyr::geom