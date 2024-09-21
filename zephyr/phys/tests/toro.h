#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/tests/test_1D.h>
#include <zephyr/phys/tests/test_2D.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;
using zephyr::geom::generator::Rectangle;

/// @class Набор тестов на распад разрыва из монографии Торо (глава 10) и (4.3.3 Numerical Tests)
/// E.F. Toro. Riemann Solvers and Numerical Methods for Fluid Dynamics.
class ToroTest : public Test1D {
public:
    IdealGas::Ptr eos_L; ///< Используемый УрС
    IdealGas::Ptr eos_R; ///< Используемый УрС

    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия


    /// @brief Конструктор
    /// @param num Номер теста 1..7
    explicit ToroTest(int num, bool multimat = false);

    /// @brief Симметрично отразить начальные условия
    void inverse();

    std::string name() const final { return "ToroTest";}

    /// @brief Левая граница области
    double xmin() const final { return 0.0; }

    /// @brief Правая граница области
    double xmax() const final { return 1.0; }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    ///@brief Получить положение разрыва
    double get_x_jump() const final { return x_jump; }


    /// @brief Начальная плотность
    double density(const Vector3d& r) const final {
        return r.x() < x_jump ? rL : rR;
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d& r) const final {
        return { r.x() < x_jump ? uL : uR, 0.0, 0.0};
    }

    /// @brief Начальное давление
    double pressure(const Vector3d& r) const final {
        return r.x() < x_jump ? pL : pR;
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d& r) const final {
        return r.x() < x_jump ? eL : eR;
    }

    /// @brief Уравнение состояния
    Eos::Ptr get_eos(const Vector3d& r) const final {
        return r.x() < x_jump ? eos_L : eos_R;
    }
};


/// @class Тесты Торо, повернутые на угол
/// Распад разрыва под углом.
class ToroTest2D {
public:
    ToroTest  test1D;
    double    angle;
    Rectangle generator;

    /// @brief Конструктор
    /// @param num Номер теста 1..7
    explicit ToroTest2D(int num, double angle, bool multimat = false)
        : test1D(num, multimat), angle(angle) {
        generator = Rectangle(-5.0, 5.0, -5.0, 5.0, false);
        generator.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
                                  .bottom=Boundary::ZOE, .top=Boundary::ZOE});
    }

    std::string name() const { return "2D Toro Test"; }

    /// @brief Конечный момент времени
    double max_time() const { return 5.0 * test1D.finish; }

    Vector3d rotate(const Vector3d& r) const {
        Vector3d res = {
                +std::cos(angle) * r.x() + std::sin(angle) * r.y(),
                -std::sin(angle) * r.x() + std::cos(angle) * r.y(),
                0.0
        };
        return res;
    }

    /// @brief Начальная плотность
    double density(const Vector3d &r) const {
        return test1D.density(rotate(r));
    }

    /// @brief Начальное давление
    double pressure(const Vector3d &r) const {
        return test1D.pressure(rotate(r));
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &r) const {
        return test1D.velocity(rotate(r));
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const {
        return test1D.energy(rotate(r));
    }

    Eos::Ptr get_eos(const Vector3d& r) const {
        return test1D.get_eos(r);
    }

    ~ToroTest2D() = default;
};

} // namespace zephyr::phys
