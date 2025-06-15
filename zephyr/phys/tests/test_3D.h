#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/geom/boundary.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/tests/ivp.h>
#include <zephyr/phys/tests/sedov.h>

namespace zephyr::phys {

using geom::generator::Cuboid;

/// @brief Абстрактный класс для трёхмерных тестов.
/// Только в прямоугольнике?
class Test3D : public IVP {
public:
    // Границы области

    virtual double xmin() const = 0;

    virtual double xmax() const = 0;

    virtual double ymin() const = 0;

    virtual double ymax() const = 0;

    virtual double zmin() const = 0;

    virtual double zmax() const = 0;

    /// @brief Граничные условия для прямоугольника
    virtual Cuboid::Boundaries boundaries() const {
        return {.left=Boundary::ZOE, .right=Boundary::ZOE,
                .bottom=Boundary::ZOE, .top=Boundary::ZOE,
                .back=Boundary::ZOE, .front=Boundary::ZOE};
    }

    /// @brief Получить название теста
    std::string name() const override { return "Test3D"; }
};

/// @brief Двумерный распад разрыва.
/// Набор тестов взят из статьи?
class SedovBlast3D : public Test3D {
public:
    // Аналитическое решение задачи
    Sedov3D exact;

    double init_time; ///< Начальный момент времени
    double finish;    ///< Конечный момент времени

    struct params {
        double gamma = 1.4;
        double rho0  = 1.0;
        double E     = 1.0;
    };

    /// @brief Конструктор
    explicit SedovBlast3D(params p);

    /// @brief Название теста
    std::string name() const final { return "3D Sedov Blast"; }

    // Границы области

    double xmin() const final { return 0.0; }

    double xmax() const final { return 1.0; };

    double ymin() const final { return 0.0; };

    double ymax() const final { return 1.0; };

    double zmin() const final { return 0.0; };

    double zmax() const final { return 1.0; };

    /// @brief Граничные условия для прямоугольника
    Cuboid::Boundaries boundaries() const final {
        return {.left=Boundary::WALL, .right=Boundary::ZOE,
                .bottom=Boundary::WALL, .top=Boundary::ZOE,
                .back=Boundary::WALL, .front=Boundary::ZOE};
    }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    // Начальные данные

    double density(const Vector3d &r) const final;

    Vector3d velocity(const Vector3d &r) const final;

    double pressure(const Vector3d &r) const final;

    // Точное решение

    double density_t(const Vector3d &r, double t) const final;

    Vector3d velocity_t(const Vector3d &r, double t) const final;

    double pressure_t(const Vector3d &r, double t) const final;
};

} // namespace zephyr::phys