#pragma once

#include <cmath>

namespace zephyr { namespace phys {

/// @brief Термодинамическая величина +
/// производные по плотности и энергии
struct dRdE {
    double value = NAN;
    double dR = NAN;
    double dE = NAN;

    dRdE(double value, double dR = NAN, double dE = NAN)
        : value(value), dR(dR), dE(dE) { }

    operator double() const {
        return value;
    };
};

/// @brief Термодинамическая величина +
/// производные по давлению и температуре
struct dPdT {
    double val = NAN;
    double dP = NAN;
    double dT = NAN;

    dPdT(double value, double dP = NAN, double dT = NAN)
            : val(value), dP(dP), dT(dT) { }

    operator double() const {
        return val;
    };
};

/// @brief Дополнительные аргументы УрС
struct Options {
    bool deriv  = false; ///< Вычислять производные
    double rho0 = NAN;   ///< Начальное приближение плотности
    double P0   = NAN;   ///< Начальное приближение давления
    double T0   = NAN;   ///< Начальное приближение температуры
};

}
}