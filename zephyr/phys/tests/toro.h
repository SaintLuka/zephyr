#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/tests/classic_test.h>

namespace zephyr {
namespace phys {

using zephyr::geom::Vector3d;

/// @class Набор тестов на распад разрыва из монографии Торо (глава 10) и (4.3.3 Numerical Tests)
/// E.F. Toro. Riemann Solvers and Numerical Methods for Fluid Dynamics.
class ToroTest : public ClassicTest {
public:

    IdealGas eos;   ///< Используемый УрС
    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double vL = 0.0, vR = 0.0;
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия


    /// @brief Конструктор
    /// @param num Номер теста 1..7
    explicit ToroTest(int num);

    /// @brief Симметрично отразить начальные условия
    void inverse();

    std::string get_name() const { return "ToroTest";}

    /// @brief Левая граница области
    double xmin() const { return 0.0; }

    /// @brief Правая граница области
    double xmax() const { return 1.0; }

    /// @brief Конечный момент времени
    double max_time() const { return finish; }

    ///@brief Получить используемый УрС
    const Eos& get_eos() const { return eos; }

    ///@brief Получить положение разрыва
    double get_x_jump() const { return x_jump; }

    /// @brief Начальная плотность
    double density(const double &x) const { return x < x_jump ? rL : rR; }

    /// @brief Начальная скорость
    Vector3d velocity(const double &x) const { return (x < x_jump) ? Vector3d(uL, vL, 0.0) : Vector3d(uR, vR, 0.0); }

    /// @brief Начальное давление
    double pressure(const double &x) const { return x < x_jump ? pL : pR; }

    /// @brief Начальная внутренняя энергия
    double energy(const double &x) const { return x < x_jump ? eL : eR; }


    /// @brief Начальная плотность
    double density(const Vector3d &r) const { return density(r.x()); }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &r) const { return velocity(r.x()); }

    /// @brief Начальное давление
    double pressure(const Vector3d &r) const { return pressure(r.x()); }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const { return energy(r.x()); }

    ~ToroTest() = default;
};

}
}
