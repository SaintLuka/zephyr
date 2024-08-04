#pragma once

#include <cmath>
#include <zephyr/phys/fractions.h>

namespace zephyr::phys {

/// @brief Термодинамическая величина +
/// производные по плотности и энергии
struct dRdE {
    double val = NAN;
    double dR = NAN;
    double dE = NAN;

    dRdE(double value, double dR = NAN, double dE = NAN)
        : val(value), dR(dR), dE(dE) { }

    operator double() const {
        return val;
    };
};

/// @brief Термодинамическая величина +
/// производные по плотности и температуре
struct dRdT {
    double val = NAN;
    double dR = NAN;
    double dT = NAN;

    dRdT(double value, double dR = NAN, double dT = NAN)
            : val(value), dR(dR), dT(dT) { }

    operator double() const {
        return val;
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

/// @brief Пара равновесных термодинамических величин
struct PairPT {
    double P;
    double T;
};

/// @brief Массив для хранения удельных объемов
using SpecVolumes = std::array<double, Fractions::max_size>;

/// @brief Дополнительные аргументы УрС
// TODO: Использовать начальные приближения для объемов
struct Options {
    bool deriv  = false; ///< Вычислять производные
    double rho0 = NAN;   ///< Начальное приближение плотности
    double P0   = NAN;   ///< Начальное приближение давления
    double T0   = NAN;   ///< Начальное приближение температуры

    /// @brief Начальные приближения для удельных объемов
    double* vols = nullptr;
};

} // namespace zephyr::phys