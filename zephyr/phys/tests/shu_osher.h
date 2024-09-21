#pragma once

#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/tests/test_1D.h>

namespace zephyr::phys {

/// @brief Тест Шу-Ошера 
/// C.-W. Shu and S. Osher. Efficient Implementation of Essentially 
/// Non-oscillatory Shock-Capturing Schemes, II (1988)
class ShuOsherTest : public Test1D {
public:
    IdealGas::Ptr eos;

    double epsilon;
    double rL, uL, pL;
    double rR, uR, pR;

    const double x_jump = -4.0;
    double finish = 1.8;

    Rectangle generator;


    /// @brief Начальные данные
    ShuOsherTest() {
        eos = IdealGas::create(1.4);

        rL = 27.0 / 7.0;
        uL = 4.0 * std::sqrt(35.0) / 9.0;
        pL = 31.0 / 3.0;

        rR  = 1.0;
        uR = 0.0;
        pR = 1.0;

        epsilon = 0.2;

        generator = Rectangle(xmin(), xmax(), -0.5, 0.5);
        generator.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
                                  .bottom=Boundary::WALL, .top=Boundary::WALL});
    }

    /// @brief Получить название теста
    std::string name() const final { return "Shu Osher"; };

    /// @brief Левая граница области
    double xmin() const final { return -5.0; }

    /// @brief Правая граница области
    double xmax() const final { return +5.0; }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    /// @brief Получить положение разрыва
    double get_x_jump() const final { return x_jump; };


    /// @brief Начальная плотность
    double density(const Vector3d& r) const final {
        return r.x() < x_jump ? rL : rR + epsilon * std::sin(5.0 * r.x());
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d& r) const final {
        return {r.x() < x_jump ? uL : uR, 0.0, 0.0};
    }

    /// @brief Начальное давление
    double pressure(const Vector3d& r) const final {
        return r.x() < x_jump ? pL : pR;;
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d& r) const final {
        return eos->energy_rp(density(r), pressure(r));
    }

    /// @brief Уравнение состояния
    Eos::Ptr get_eos(const Vector3d& r) const final {
        return eos;
    }
};

} // namespace zephyr::phys
