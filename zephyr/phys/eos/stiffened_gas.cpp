#include <string>
#include <stdexcept>
#include <cmath>

#include <zephyr/phys/eos/stiffened_gas.h>

namespace zephyr::phys {

StiffenedGas::StiffenedGas(double gamma, double P0, double e0, double Cv, double T0)
    : gamma(gamma), P0(P0), e0(e0), Cv(Cv), T0(T0) { }

StiffenedGas::StiffenedGas(const std::string &name) {
    // @formatter:off
    P0 = 0.0;
    e0 = 0.0;
    Cv = 1.0;
    T0 = 0.0;

    if (name == "Air") {
        gamma = 1.4;
        Cv    = 718.0_J_kgK;
    }
    else if (name == "SF6") {
        gamma = 1.1074;
        Cv    = 567.0_J_kgK;
    }
    else if (name == "Water") {
        gamma = 4.4;
        Cv    = 4160.0_J_kgK;
        P0    = 600.0_MPa;
        e0    = 0.0;
        T0    = 230.7; // 10^5 Па, 1000 кг/м^3, 0 °C
    }
    else if (name == "Water2") {
        gamma = 3.0;
        Cv    = 4160.0_J_kgK;
        P0    = 853.3_MPa;
        e0    = -1.148_MJ_kg;
    }
    else if (name == "Copper" || name == "Cu") {
        gamma = 4.0;
        Cv    = 384.0_J_kgK;
        P0    = 341.0;
        e0    = 0.0;
    }
    else if (name == "Lead" || name == "Pb") {
        // Из параметров для Mie-Gruneisen
        // при (P = 0.0, rho = rho0)
        gamma = 3.77;
        Cv    = 121.0_J_kgK;
        P0    = 11.8_GPa;
        e0    = -1.42_MJ_kg;
        T0    = -2832.57;
    }
    else if (name == "Gas") {
        // Тестовый газ (аналог Air)
        gamma = 1.4;
        Cv    = 2.5;
    }
    else if (name == "HeavyGas") {
        // Тестовый тяжелый газ (аналог SF6)
        gamma = 1.1;
        Cv    = 2.5;
    }
    else if (name == "Liquid") {
        // Тестовая жидкость (аналог воды)
        gamma = 4.4;
        Cv    = 20.0;
        P0    = 6.8e4;
    }
    else if (name == "Solid") {
        // Тестовый металл (аналог свинца)
        gamma = 3.8;
        Cv    = 0.5;
        P0    = 1.4e5;
        e0    = -15.0;
        T0    = -7.0;
    }
    else {
        throw std::runtime_error("Unknown stiffened gas '" + std::string(name) + "'");
    }
    // @formatter:on
}

dRdE StiffenedGas::pressure_re(double rho, double eps, const Options& options) const {
    dRdE res {(gamma - 1.0) * rho * (eps - e0) - gamma * P0 };
    if (options.deriv) {
        res.dR = (gamma - 1.0) * (eps - e0);
        res.dE = (gamma - 1.0) * rho;
    }
    return res;
}

dRdT StiffenedGas::pressure_rT(double rho, double T, const Options& options) const {
    dRdT P = (gamma - 1.0) * Cv * rho * (T - T0) - P0;
    if (options.deriv) {
        P.dR = (gamma - 1.0) * Cv * (T - T0);
        P.dT = (gamma - 1.0) * Cv * rho;
    }
    return P;
}

dRdT StiffenedGas::energy_rT(double rho, double T, const Options& options) const {
    dRdT e = e0 + Cv * (T - T0) + P0 / rho;
    if (options.deriv) {
        e.dR = -P0 / (rho * rho);
        e.dT = Cv;
    }
    return e;
}

double StiffenedGas::energy_rP(double rho, double P, const Options& options) const {
    return e0 + (P + gamma * P0) / ((gamma - 1.0) * rho);
}

double StiffenedGas::sound_speed_re(double rho, double eps, const Options& options) const {
    return std::sqrt(gamma * (gamma - 1.0) * ((eps - e0) - P0 / rho));
}

double StiffenedGas::sound_speed_rP(double rho, double P, const Options& options) const {
    return std::sqrt(gamma * (P + P0) / rho);
}

double StiffenedGas::temperature_rP(double rho, double P, const Options& options) const {
    return T0 + (P + P0) / ((gamma - 1.0) * Cv * rho);
}

dPdT StiffenedGas::volume_PT(double P, double T, const Options& options) const {
    dPdT res{((gamma - 1.0) * Cv * (T - T0)) / (P + P0)};
    if (options.deriv) {
        res.dP = -res.val / (P + P0);
        res.dT = +res.val / (T - T0);
    }
    return res;
}

inline double sqr(double x) { return x * x; }

dPdT StiffenedGas::energy_PT(double P, double T, const Options& options) const {
    dPdT res{e0 + ((P + gamma * P0) / (P + P0)) * Cv * (T - T0)};
    if (options.deriv) {
        res.dP = (1.0 - gamma) * P0 * Cv * (T - T0) / sqr(P + P0);
        res.dT = ((P + gamma * P0) / (P + P0)) * Cv;
    }
    return res;
}

StiffenedGas StiffenedGas::stiffened_gas(double rho, double P, const Options& options) const {
    return StiffenedGas(gamma, P0, e0, Cv, T0);
}

double StiffenedGas::min_pressure() const {
    return -P0;
}

} // namespace zephyr::phys