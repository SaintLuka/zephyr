#pragma once

#include <zephyr/geom/curves/interpolant.h>

namespace zephyr::geom::curves {

/// @brief Интерполяционный линейный сплайн y(x).
class LinearSpline : public Interpolant {
public:
    /// @brief Конструктор линейного сплайна y(x)
    /// @param x, y Аргумент/значение в узлах
    /// @param left, right Тип экстрапляции на границах
    LinearSpline(const std::vector<double>& x,
            const std::vector<double>& y,
            SplineBound left  = SplineBound::Free,
            SplineBound right = SplineBound::Free);

    /// @brief Основная функция. Значение функции от аргумента
    double get(double x) const override;
};


/// @brief Ломаная линия на наборе точек.
/// Параметрическая кривая x(t), y(t), z(t), параметр t ∈ [0, 1].
class PLinearSpline : public PInterpolant {
public:
    /// @brief Конструктор
    /// @param xs, ys Узлы ломаной
    /// @param left, right Тип экстрапляции на границах
    PLinearSpline(const std::vector<double> &xs,
                  const std::vector<double> &ys,
                  SplineBound left = SplineBound::Crop,
                  SplineBound right = SplineBound::Crop);

    /// @brief Конструктор
    /// @param xs, ys, zs Узлы ломаной
    /// @param left, right Тип экстрапляции на границах
    PLinearSpline(const std::vector<double> &xs,
                  const std::vector<double> &ys,
                  const std::vector<double> &zs,
                  SplineBound left = SplineBound::Crop,
                  SplineBound right = SplineBound::Crop);

    /// @brief Конструктор ломаной линии.
    /// @param vs Точки, на которых строится ломаная
    PLinearSpline(const std::vector<Vector3d> &vs,
                  SplineBound left = SplineBound::Crop,
                  SplineBound right = SplineBound::Crop);


    /// @brief Получить точку на ломаной
    double x(double t) const override;

    /// @brief Получить точку на ломаной
    double y(double t) const override;

    /// @brief Получить точку на ломаной
    double z(double t) const override;

private:
    void build(const std::vector<double> &xs,
               const std::vector<double> &ys,
               const std::vector<double> &zs,
               SplineBound left, SplineBound right);
};

} // namespace zephyr::geom::curves