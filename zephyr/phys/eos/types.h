#pragma once

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
/// производные по плотности и давлению
struct dRdP {
    double value = NAN;
    double dR = NAN;
    double dP = NAN;

    dRdP(double value, double dR = NAN, double dP = NAN)
            : value(value), dR(dR), dP(dP) { }

    operator double() const {
        return value;
    };
};

/// @brief Термодинамическая величина +
/// производные по давлению и температуре
struct dPdT {
    double value = NAN;
    double dP = NAN;
    double dT = NAN;

    dPdT(double value, double dP = NAN, double dT = NAN)
            : value(value), dP(dP), dT(dT) { }

    operator double() const {
        return value;
    };
};

/// @brief Дополнительные аргументы УрС
struct Options {
    bool deriv  = false; ///< Вычислять производные
    double rho0 = NAN;   ///< Начальное приближение плотности
};

}
}