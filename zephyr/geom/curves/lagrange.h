#pragma once

#include <zephyr/geom/curves/interpolant.h>

namespace zephyr::geom::curves {

/// @brief Интерполяционный полином Лагранжа y(x)
class Lagrange : public Interpolant {
public:
    /// @brief Конструктор полинома Лагранжа y(x)
    /// @param x, y Аргумент/значение в узлах
    /// @param left, right Тип экстрапляции на границах
    Lagrange(
            const std::vector<double>& x,
            const std::vector<double>& y,
            SplineBound left  = SplineBound::Free,
            SplineBound right = SplineBound::Free);

    /// @brief Основная функция. Значение функции от аргумента
    double get(double x) const override;
};


/// @brief Полином Лагранжа по каждой координате.
/// Параметрическая кривая x(t), y(t), z(t), параметр t ∈ [0, 1].
class PLagrange : public PInterpolant {
public:
    /// @brief Конструктор
    /// @param xs, ys Узлы сплайна
    /// @param left, right Тип экстрапляции на границах
    PLagrange(const std::vector<double> &xs,
              const std::vector<double> &ys,
              SplineBound left = SplineBound::Crop,
              SplineBound right = SplineBound::Crop,
              Parametrization param = Parametrization::Chebyshev);

    /// @brief Конструктор
    /// @param xs, ys, zs Узлы сплайна
    /// @param left, right Тип экстрапляции на границах
    PLagrange(const std::vector<double> &xs,
              const std::vector<double> &ys,
              const std::vector<double> &zs,
              SplineBound left = SplineBound::Crop,
              SplineBound right = SplineBound::Crop,
              Parametrization param = Parametrization::Chebyshev);

    /// @brief Конструктор ломаной линии.
    /// @param vs Точки, на которых строится сплайн
    PLagrange(const std::vector<Vector3d> &vs,
              SplineBound left = SplineBound::Crop,
              SplineBound right = SplineBound::Crop,
              Parametrization param = Parametrization::Chebyshev);


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
};

} // namespace zephyr::geom::curves