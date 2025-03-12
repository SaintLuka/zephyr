#pragma once

#include <zephyr/geom/vector.h>

namespace zephyr::geom {

/// @namespace zephyr::geom::obj
/// @brief Простейшие геометрические объекты: отрезок, линия, луч, плоскость
namespace obj {

/// @brief Отрезок
struct segment {
    Vector3d v1;  ///< Начальная точка
    Vector3d v2;  ///< Конечная точка

    /// Параметризация: r(t) = v1 + (v2 - v1) * t
    Vector3d get(double t) const {
        return v1 + (v2 - v1) * t;
    }

    /// @brief Направляющая
    Vector3d tau() const {
        return v2 - v1;
    }

    /// @brief Длина
    double length() const {
        return (v2 - v1).norm();
    }
};

/// @brief Прямая
struct line {
    Vector3d p;    ///< Точка прямой
    Vector3d tau;  ///< Направляющая

    /// @brief Параметризация: r(t) = p + tau * t
    Vector3d get(double t) const {
        return p + tau * t;
    }
};

/// @brief Луч
struct ray {
    Vector3d p;    ///< Начало луча
    Vector3d tau;  ///< Направление луча

    /// @brief Параметризация r(t) = p + tau * t
    Vector3d get(double t) const {
        return p + tau * t;
    }
};

/// @brief Плоскость
/// @details Множество точек r: (r - p, n) = 0.
///          Считается, что нормаль n -- внешняя.
struct plane {
    Vector3d p;   ///< Точка плоскости
    Vector3d n;   ///< Внешняя нормаль плоскости

    /// @brief Параметризация в двумерном случае (прямая)
    /// Нормаль n - правая при движении вдоль прямой.
    Vector3d get(double t) const {
        return {p.x() - n.y() * t, p.y() + n.x() * t, 0.0};
    }

    /// @brief Под плоскостью? (n -- внешняя нормаль)
    bool under(const Vector3d& v) {
        return (v - p).dot(n) < 0.0;
    }

    /// @brief Положение точки относительно плоскости
    /// @return -1: под, 0: на, +1: над
    int position(const Vector3d& v) {
        double val = (v - p).dot(n);
        return val == 0.0 ? 0 : (val < 0.0 ? -1 : +1);
    }
};

/// @brief Окружность
struct circle {
    Vector3d c;  ///< Центр окружности
    double r;    ///< Радиус окружности
};

} // namespace obj


/// @namespace zephyr::geom::intersection2D
/// @brief Пересечения на плоскости
/// @details В функциях данного пространства имен подразумевается,
/// что объекты являются двумерными, то есть координата z = 0.
/// В данном случае структура plane подразумевается прямой.
namespace intersection2D {

/// @brief Существует пересечение луча и отрезка?
/// @details Пересечение со второй точкой отрезка не учитывается
bool exist(const obj::ray& ray, const obj::segment& seg);

/// @brief Существует пересечение плоскости и отрезка?
bool exist(const obj::plane& plane, const obj::segment& seg);

/// @brief Найти пересечение двух линий без всяких вопросов
/// @details Для параллельных линий, вероятнее всего, выдаст NAN
Vector3d find_fast(const obj::line& line1, const obj::line& line2);

/// @brief Найти пересечение двух линий без всяких вопросов
/// @details Для параллельных линий, вероятнее всего, выдаст NAN
/// Точка пересечения может лежать вне отрезков
Vector3d find_fast(const obj::segment& seg1, const obj::segment& seg2);

/// @brief Найти пересечение двух линий без всяких вопросов
/// @details Для параллельных линий, вероятнее всего, выдаст NAN
/// Точка пересечения может лежать вне отрезка
Vector3d find_fast(const obj::line& line, const obj::segment& seg);

/// @brief Найти пересечение плоскости и отрезка без всяких вопросов
/// @details Для параллельных, вероятнее всего, выдаст NAN
/// Точка пересечения может лежать вне отрезка
Vector3d find_fast(const obj::plane& plane, const obj::line& line);

/// @brief Найти пересечение плоскости и отрезка без всяких вопросов
/// @details Для параллельных, вероятнее всего, выдаст NAN
/// Точка пересечения может лежать вне отрезка
Vector3d find_fast(const obj::plane& plane, const obj::segment& seg);


/// @brief Пересечение окружности и отрезка
/// @details Пусть прямая задана параметрически: r(t) = p + tau * t
/// Если пересечения с окружностью существуют, тогда находятся
/// параметры t1, t2 и соответствующие точки пересечения p1, p2
/// Пересечение по касательной игнорируется.
struct circle_segment_intersection {
    bool exist;       ///< Пересечения существуют?
    double t1, t2;    ///< Параметризация прямой t1 < t2
    Vector3d p1, p2;  ///< Точки пересечения
};

/// @brief Пересечение окружности и прямой
circle_segment_intersection find(const obj::circle& circle, const obj::line& line);

/// @brief Пересечение окружности и прямой (пересечения вне отрезка считаются)
circle_segment_intersection find(const obj::circle& circle, const obj::segment& seg);


} // namespace intersection2D



/// @namespace zephyr::geom::intersection3D
/// @brief Пересечения в пространстве
namespace intersection3D {

} // namespace intersection3D

} // namespace zephyr::geom