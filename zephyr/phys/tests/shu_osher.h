#pragma once

#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/tests/test_1D.h>

namespace zephyr::phys {

/// @brief Тест Шу-Ошера 
/// C.-W. Shu and S. Osher. Efficient Implementation of Essentially 
/// Non-oscillatory Shock-Capturing Schemes, II (1988)
class ShuOsherTest : public Test1D {
public:

    /// @brief Начальные данные
    ShuOsherTest() {
        m_eos = IdealGas::create(1.4);

        rL  = 27.0 / 7.0;
        uL = 4.0 * std::sqrt(35.0) / 9.0;
        pL = 31.0 / 3.0;

        rR  = 1.0;
        uR = 0.0;
        pR = 1.0;

        epsilon = 0.2;
    }

    /// @brief Получить название теста
    std::string get_name() const { return "Shu Osher"; };

    /// @brief Используемое уравнение состояния
    Eos::Ptr get_eos() const final { return m_eos; }

    /// @brief Левая граница области
    double xmin() const final { return -5.0; }

    /// @brief Правая граница области
    double xmax() const final { return +5.0; }

    /// @brief Конечный момент времени
    double max_time() const final { return 2.0; }

    /// @brief Получить положение разрыва
    double get_x_jump() const final { return x_jump; };


    /// @brief Начальная плотность
    double density(double x) const final {
        return x < x_jump ? rL : rR + epsilon * std::sin(5.0 * x);
    }

    /// @brief Начальная скорость
    Vector3d velocity(double x) const final {
        return {x < x_jump ? uL : uR, 0.0, 0.0};
    }

    /// @brief Начальное давление
    double pressure(double x) const final {
        return x < x_jump ? pL : pR;;
    }

    /// @brief Начальная внутренняя энергия
    double energy(double x) const final {
        return m_eos->energy_rp(density(x), pressure(x));
    }


    /// @brief Начальная плотность
    double density(const Vector3d& r) const final {
        return density(r.x());
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d& r) const final {
        return velocity(r.x());
    }

    /// @brief Начальное давление
    double pressure(const Vector3d& r) const final {
        return pressure(r.x());
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d& r) const final {
        return energy(r.x());
    }

private:
    IdealGas::Ptr m_eos;

    double epsilon;
    double rL, uL, pL;
    double rR, uR, pR;

    const double x_jump = -4.0;
};

} // namespace zephyr::phys
