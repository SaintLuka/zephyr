#pragma once

#include <zephyr/geom/geom.h>

namespace zephyr { namespace geom {

/// @class Задание квадратичного сплайна в параметрическом виде.
/// Отображение параметра -1 < xi < 1 на кривую, проходящую
/// через три точки v1, vc, v2
struct Mapping1D {

    Vector3d v1; ///< Начальная точка (соответствует t = -1)
    Vector3d vc; ///< Центральная точка (при t = 0)
    Vector3d v2; ///< Конечная точка (соответствует t = +1)

    /// @brief Конструктор простого отрезка
    explicit Mapping1D(const ShortList1D &vs);

    /// @brief Конструктор квадратичного сплайна
    explicit Mapping1D(const LargeList1D &vs);

    /// @brief Получить точку на сплайне без создания экземпляра класса
    static Vector3d get(const Vector3d &v1, const Vector3d &vc,
                      const Vector3d &v2, double xi);

    /// @brief Получить касательный вектор без создания экземпляра класса
    static Vector3d tangent(const Vector3d &v1, const Vector3d &vc,
                          const Vector3d &v2, double xi);

    /// @brief Получить точку на сплайне, соответствующую параметру xi
    Vector3d operator()(double xi) const;

    /// @brief Произвольная нормаль, в двумерных задачах является правой
    /// нормалью к направленному сплайну
    Vector3d normal(double xi) const;

    /// @brief Якобиан отображения. Фактически модуль производной.
    double Jacobian(double xi) const;
};


/// @class Квадратичное отображение квадрата -1 < xi < 1, -1 < eta < 1
/// на произвольный четырехугольник.
struct Mapping2D {

    Mapping1D SB; ///< Нижний сплайн
    Mapping1D SH; ///< Горизонтальный сплайн
    Mapping1D ST; ///< Верхний сплайн

    /// @brief Создание линейного отображения на четырехугольник
    explicit Mapping2D(const ShortList2D& vs);

    /// @brief Создание квадратичного отображения на четырехугольник
    explicit Mapping2D(const LargeList2D& vs);

    /// @brief Получить точку четырехугольника по значениям параметров
    Vector3d operator()(double xi, double eta) const;

    /// @brief Нормаль
    Vector3d normal(double xi, double eta) const;

    /// @brief Якобиан отображения
    double Jacobian(double xi, double eta) const;
};


/// @class Функтор, предоставляющий отображение куба
/// |xi| < 1, |eta| < 1, |chi| < 1 на шестигранник с восемью вершинами
/// (топологический куб)
struct Mapping3D {

    /// @brief Набор вершин ячейки
    const ShortList3D& vs;

    /// @brief Конструктор простой ячейки
    explicit Mapping3D(const ShortList3D& vs);

    /// @brief Непосредственно отображение
    Vector3d operator()(double xi, double eta, double chi);

    /// @brief Якобиан отображения
    double Jacobian(double xi, double eta, double chi);
};

} // geom
} // zephyr