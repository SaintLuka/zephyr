#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr { namespace phys {

/// @class Классический одномерный тест Сода
class SodTest {
public:
    IdealGas eos;   ///< Используемый УрС
    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия

    /// @brief Конструктор
    SodTest() : eos(1.4) {
        rL = 1.0; rR = 0.125;
        pL = 1.0; pR = 0.1;
        uL = 0.0; uR = 0.0;
        eL = eos.energy_rp(rL, pL);
        eR = eos.energy_rp(rR, pR);

        x_jump = 0.5;
        finish = 0.2;
    }

    /// @brief Симметрично отразить начальные условия
    void inverse() {
        std::swap(rL, rR);
        std::swap(uL, uR);
        std::swap(pL, pR);
        std::swap(eL, eR);
        uL *= -1.0;
        uR *= -1.0;
    }


    /// @brief Левая граница области
    double xmin() const { return 0.0; }

    /// @brief Правая граница области
    double xmax() const { return 1.0; }

    /// @brief Конечный момент времени
    double max_time() const { return finish; }


    /// @brief Начальная плотность
    double density(const double &x) const { return x < x_jump ? rL : rR; }

    /// @brief Начальная скорость
    Vector3d velocity(const double &x) const { return {x < x_jump ? uL : uR, 0.0, 0.0}; }

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

};

}
}
