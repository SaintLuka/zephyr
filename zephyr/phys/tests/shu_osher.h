#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr { namespace phys {

/// @brief Тест Шу-Ошера 
/// C.-W. Shu and S. Osher. Efficient Implementation of Essentially 
/// Non-oscillatory Shock-Capturing Schemes, II (1988)
class ShuOsherTest {
public:

    /// @brief Начальные данные
    ShuOsherTest() : m_eos(1.4) {
        rL  = 27.0 / 7.0;
        uL = 4.0 * std::sqrt(35.0) / 9.0;
        pL = 31.0 / 3.0;

        rR  = 1.0;
        uR = 0.0;
        pR = 1.0;

        epsilon = 0.2;
    }

    /// @brief Используемое уравнение состояния
    IdealGas& eos() { return m_eos; }

    /// @brief Левая граница области
    double xmin() const { return -5.0; }

    /// @brief Правая граница области
    double xmax() const { return +5.0; }

    /// @brief Конечный момент времени
    double max_time() const { return 2.0; }


    /// @brief Начальная плотность
    double density(const Vector3d &r) const {
        return r.x() < x_jump ? rL : rR + epsilon * std::sin(5.0 * r.x());
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &r) const {
        return {r.x() < x_jump ? uL : uR, 0.0, 0.0};
    }

    /// @brief Начальное давление
    double pressure(const Vector3d &r) const {
        return r.x() < x_jump ? pL : pR;;
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const {
        return m_eos.energy_rp(density(r), pressure(r));
    }

private:
    IdealGas m_eos;

    double epsilon;
    double rL, uL, pL;
    double rR, uR, pR;

    const double x_jump = -4.0;
};

}
}
