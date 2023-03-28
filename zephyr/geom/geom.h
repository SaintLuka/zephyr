#pragma once

#include <array>
#include <vector>

#include <zephyr/geom/vector.h>


namespace zephyr { namespace geom {

/// @brief Представление отрезка
using ShortList1D = std::array<Vector3d, 2>;

/// @brief Точки квадратичного сплайна
using LargeList1D = std::array<Vector3d, 3>;

/// @brief Представление четырехугольника
using ShortList2D = std::array<Vector3d, 4>;

/// @brief Квадратичное отображение на четырехугольник
using LargeList2D = std::array<Vector3d, 9>;

/// @brief Представление шестигранника
using ShortList3D = std::array<Vector3d, 8>;

/// @brief Квадратичное отображение на шестигранник (not used)
using LargeList3D = std::array<Vector3d, 27>;

/// @brief Список вершин произвольной длины
/// для описания многоугольников
using VerticesList = std::vector<Vector3d>;


/// @brief Длина простой одномерной грани, равна длине отрезка,
/// соединяющего две точки.
double length(const ShortList1D& vs);

/// @brief "Длина" криволинейной одномерной грани. Равна | int vec(n) dl |,
/// то есть модулю криволинейного интеграла 2-го рода.
double length(const LargeList1D& vs);

/// @brief Внешняя нормаль к простой одномерной грани.
/// Вершины должны располагаться в плоскости (x, y).
/// @param center Центр ячейки
Vector3d normal(const ShortList1D& vs, const Vector3d& center);

/// @brief Внешняя нормаль к криволинейной одномерной грани.
/// Вершины должны располагаться в плоскости (x, y).
/// @param center Центр ячейки
Vector3d normal(const LargeList1D& vs, const Vector3d& center);

/// @brief Внешняя нормаль к простой двумерной грани
/// @param center Центр ячейки
Vector3d normal(const ShortList2D& vs, const Vector3d& center);

/// @brief "Площадь" простой двумерной грани. Равна | int vec(n) dS |,
/// то есть модулю поверхностного интеграла 2-го рода.
double area(const ShortList2D &vs);

/// @brief Площадь криволинейной двумерной ячейки.
/// Вершины должны располагаться в плоскости (x, y).
double area(const LargeList2D &vs);

/// @brief Площадь произвольного многоугольника
/// Вершины должны располагаться в плоскости (x, y).
double area(const VerticesList& vs);

/// @brief Объем обычной двумерной ячейки в осесимметричной постновке
/// Вершины должны располагаться в плоскости (x, y), где y - радиус
double volume_as(const ShortList2D &vs);

/// @brief Объем криволинейной двумерной ячейки в осесимметричной постновке
/// Вершины должны располагаться в плоскости (x, y), где y - радиус
double volume_as(const LargeList2D &vs);

/// @brief Объем простой трехмерной ячейки
double volume(const ShortList3D &vs);

/// @brief Центр отрезка
Vector3d center(const ShortList1D &vs);

/// @brief Центр отрезка
Vector3d center(const LargeList1D &vs);

/// @brief Центр простой двумерной ячейки
Vector3d center(const ShortList2D &vs);

/// @brief Центр криволинейной двумерной ячейки
Vector3d center(const LargeList2D &vs);

/// @brief Центр произвольного многоугольника (среднее вершин)
Vector3d center(const VerticesList& vs);

/// @brief Центр трехмерной ячейки
Vector3d center(const ShortList3D &vs);

/// @brief Барицентр простой двумерной ячейки.
/// Вершины должны располагаться в плоскости (x, y).
/// @param area Площадь ячейки
Vector3d centroid(const ShortList2D &vs, double area = 0.0);

/// @brief Барицентр произвольного многоугольника.
/// Вершины должны располагаться в плоскости (x, y).
/// @param area Площадь многоугольника
Vector3d centroid(const VerticesList &vs, double area = 0.0);

/// @brief Барицентр криволинейной двумерной ячейки.
/// Вершины должны располагаться в плоскости (x, y).
/// @param area Площадь криволинейной ячейки
Vector3d centroid(const LargeList2D &vs, double area = 0.0);

/// @brief Барицентр трехмерной ячейки
/// @param volume Объем ячейки
Vector3d centroid(const ShortList3D &vs, double volume = 0.0);

/// @brief Проверка всех вышеперечисленных функций
/// путем сравнения результатов численного интегрирования
void checkout();

} // geom
} // zephyr