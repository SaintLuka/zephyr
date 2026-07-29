#pragma once
#include <tuple>
#include <memory>

#include <zephyr/geom/vector.h>

namespace zephyr::geom {

/// @brief Топография, уровень дна
class IBed {
public:
    using Ptr = std::shared_ptr<IBed>;
    using Ref = const std::shared_ptr<IBed>&;

    virtual ~IBed() = default;

    /// @brief Высота дна и градиент в зависимости от координаты
    virtual std::tuple<double, Vector2d> get(const Vector3d& v) const {
        throw std::runtime_error("IBed::get() is not implemented");
    }
};

/// @brief Ровное дно
class ConstBed : public IBed {
    std::tuple<double, Vector2d> res;

public:
    explicit ConstBed(double bottom)
        : res({bottom, Vector2d::Zero()}) {}

    /// @param bottom Уровень дна
    static IBed::Ptr create(double bottom) {
        return std::make_shared<ConstBed>(bottom);
    }

    /// @brief Высота дна и градиент в зависимости от координаты
    std::tuple<double, Vector2d> get(const Vector3d& v) const final {
        return res;
    }
};

/// @brief Дно в виде ступеньки
class StepBed : public IBed {
    double left, right, x0;

public:
    explicit StepBed(double left, double right, double x0 = 0.0)
        : left(left), right(right), x0(x0) {}

    /// @param
    static IBed::Ptr create(double lvl1, double lvl2, double x0 = 0.0) {
        return std::make_shared<StepBed>(lvl1, lvl2, x0);
    }

    /// @brief Высота дна и градиент в зависимости от координаты
    std::tuple<double, Vector2d> get(const Vector3d& v) const final {
        return {(v.x() < x0 ? left : right), Vector2d::Zero()};
    }
};

/// @brief Плоское дно: A*x + B*y + C
class PlainBed : public IBed {
    double A, B, C;
public:
    PlainBed(double A, double B, double C)
        : A(A), B(B), C(C) { }

    /// @param A, B, C Коэффициент плоскости: z = A*x + B*y + C
    static IBed::Ptr create(double A, double B, double C) {
        return std::make_shared<PlainBed>(A, B, C);
    }

    /// @brief Высота дна и градиент в зависимости от координаты
    std::tuple<double, Vector2d> get(const Vector3d& v) const final {
        double res = A * v.x() + B * v.y() + C;
        return {res, Vector2d{A, B}};
    }
};

/// @brief Параболическое дно: coeff * (x^2 + y^2) + bottom
class ParabolicBed : public IBed {
    double coeff, bottom;
public:
    ParabolicBed(double coeff, double bottom)
        : coeff(coeff), bottom(bottom) { }

    /// @param coeff Коэффициент параболоида
    /// @param bottom Уровень нижней точки
    static IBed::Ptr create(double coeff, double bottom) {
        return std::make_shared<ParabolicBed>(coeff, bottom);
    }

    /// @brief Высота дна и градиент в зависимости от координаты
    std::tuple<double, Vector2d> get(const Vector3d& v) const final {
        double res = coeff * v.squaredNorm() + bottom;
        return {res, Vector2d{2.0*coeff*v.x(), 2.0*coeff*v.y()}};
    }
};

/// @brief Выбоина с выходом на ровное дно
class PitBed : public IBed {
    double radius, bottom, top;

public:
    PitBed(double radius, double bottom, double top)
        : radius(radius), bottom(bottom), top(top) {}

    /// @param radius Радиус выбоины
    /// @param bottom Нижняя точка
    /// @param top Уровень плоского дна
    static IBed::Ptr create(double radius, double bottom, double top) {
        return std::make_shared<PitBed>(radius, bottom, top);
    }

    /// @brief Высота дна и градиент в зависимости от координаты
    std::tuple<double, Vector2d> get(const Vector3d& v) const final {
        double r = v.norm();
        double res = top;
        Vector2d slope = Vector2d::Zero();
        if (r < radius) {
            res = 0.5 * (bottom + top) - 0.5 * (top - bottom) * std::cos(M_PI * r / radius);
            slope = {NAN, NAN};
        }
        return {res, slope};
    }
};

}