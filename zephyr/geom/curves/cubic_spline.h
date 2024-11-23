#pragma once

#include <zephyr/geom/curves/interpolant.h>

namespace zephyr::geom::curves {

/// @brief Интерполяционный кубический сплайн y(x).
class CubicSpline : public Interpolant {
public:
    /// @brief Конструктор кубического сплайна y(x)
    /// @param x, y Аргумент/значение в узлах
    /// @param left, right Тип экстрапляции на границах
    CubicSpline(const std::vector<double> &x,
                const std::vector<double> &y,
                SplineBound left = SplineBound::Free,
                SplineBound right = SplineBound::Free);

    /// @brief Основная функция. Значение функции от аргумента
    double get(double x) const override;

private:
    // Параметры сплайнов (размер на единицу меньше числа узлов)
    std::vector<double> m_b;
    std::vector<double> m_c;
    std::vector<double> m_d;
};

/// @brief Параметрический кубический сплайн.
/// Параметрическая кривая x(t), y(t), z(t), параметр t ∈ [0, 1].
class PCubicSpline : public PInterpolant {
public:
    /// @brief Конструктор
    /// @param xs, ys Узлы сплайна
    /// @param left, right Тип экстрапляции на границах
    PCubicSpline(const std::vector<double> &xs,
                 const std::vector<double> &ys,
                 SplineBound left = SplineBound::Crop,
                 SplineBound right = SplineBound::Crop,
                 Parametrization param = Parametrization::Chord);

    /// @brief Конструктор
    /// @param xs, ys, zs Узлы сплайна
    /// @param left, right Тип экстрапляции на границах
    PCubicSpline(const std::vector<double> &xs,
                 const std::vector<double> &ys,
                 const std::vector<double> &zs,
                 SplineBound left = SplineBound::Crop,
                 SplineBound right = SplineBound::Crop,
                 Parametrization param = Parametrization::Chord);

    /// @brief Конструктор ломаной линии.
    /// @param vs Точки, на которых строится сплайн
    PCubicSpline(const std::vector<Vector3d> &vs,
                 SplineBound left = SplineBound::Crop,
                 SplineBound right = SplineBound::Crop,
                 Parametrization param = Parametrization::Chord);


    /// @brief Получить точку на кривой
    double x(double t) const override;

    /// @brief Получить точку на кривой
    double y(double t) const override;

    /// @brief Получить точку на кривой
    double z(double t) const override;

private:
    void build(const std::vector<double> &xs,
               const std::vector<double> &ys,
               const std::vector<double> &zs,
               SplineBound left, SplineBound right,
               Parametrization param);

    // Параметры сплайнов (размер на единицу меньше)
    std::vector<double> m_xb;
    std::vector<double> m_yb;
    std::vector<double> m_zb;

    std::vector<double> m_xc;
    std::vector<double> m_yc;
    std::vector<double> m_zc;

    std::vector<double> m_xd;
    std::vector<double> m_yd;
    std::vector<double> m_zd;
};

} // namespace zephyr::geom::curves