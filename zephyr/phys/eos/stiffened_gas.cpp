#include <string>
#include <stdexcept>
#include <cmath>

#include <zephyr/phys/eos/stiffened_gas.h>

namespace zephyr { namespace phys {

StiffenedGas::StiffenedGas(double gamma, double p_inf, double eps_0, double Cv)
    : gamma(gamma), P0(p_inf), eps_0(eps_0), Cv(Cv) { }

StiffenedGas::StiffenedGas(const std::string &name) {
    // @formatter:off

    if (name == "Air") {
        gamma = 1.4;
        Cv    = 0.0;
        P0    = 0.0;
        eps_0 = 0.0;
    }
    else if (name == "SF6") {
        gamma = 1.076;
        Cv    = 96.6_J_kgK;
        P0    = 0.0;
        eps_0 = 0.0;
    }
    else if (name == "Water") {
        gamma = 4.4;
        Cv    = 1400.0_J_kgK;
        P0    = 600.0_MPa;
        eps_0 = 0.0;
    }
    else if (name == "Water2") {
        gamma = 3.0;
        Cv    = 1400.0_J_kgK;
        P0    = 853.3_MPa;
        eps_0 = -1.148_MJ_kg;
    }
    else if (name == "Copper") {
        gamma = 4.0;
        Cv    = 384.0_J_kgK;
        P0    = 341.0;
        eps_0 = 0.0;
    }
    else {
        throw std::runtime_error("Unknown stiffened gas '" + std::string(name) + "'");
    }
    // @formatter:on
}

dRdE StiffenedGas::pressure_re(double rho, double eps, const Options& options) const {
    dRdE res {(gamma - 1.0) * rho * (eps - eps_0) - gamma * P0 };
    if (options.deriv) {
        res.dR = (gamma - 1.0) * (eps - eps_0);
        res.dE = (gamma - 1.0) * rho;
    }
    return res;
}

double StiffenedGas::energy_rp(double rho, double P, const Options& options) const {
    return eps_0 + (P + gamma * P0) / ((gamma - 1.0) * rho);
}

double StiffenedGas::sound_speed_re(double rho, double eps, const Options& options) const {
    return std::sqrt(gamma * (gamma - 1.0) * ((eps - eps_0) - P0 / rho));
}

double StiffenedGas::sound_speed_rp(double rho, double P, const Options& options) const {
    return std::sqrt(gamma * (P + P0) / rho);
}

double StiffenedGas::pressure_rt(double rho, double T, const Options& options) const {
    return (gamma - 1.0) * Cv * rho * T - P0;
}

double StiffenedGas::temperature_rp(double rho, double P, const Options& options) const {
    return (P + P0) / ((gamma - 1.0) * Cv * rho);
}

dPdT StiffenedGas::volume_pt(double P, double T, const Options& options) const {
    dPdT res{((gamma - 1.0) * Cv * T) / (P + P0)};
    if (options.deriv) {
        res.dP = -res.val / (P + P0);
        res.dT = +res.val / T;
    }
    return res;
}

inline double sqr(double x) { return x * x; }

dPdT StiffenedGas::energy_pt(double P, double T, const Options& options) const {
    dPdT res{eps_0 + ((P + gamma * P0) / (P + P0)) * Cv * T};
    if (options.deriv) {
        res.dP = (1.0 - gamma) * P0 * Cv * T / sqr(P + P0);
        res.dT = ((P + gamma * P0) / (P + P0)) * Cv;
    }
    return res;
}

StiffenedGas StiffenedGas::stiffened_gas(double rho, double P, const Options& options) const {
    return StiffenedGas(gamma, P0, eps_0, Cv);
}

}
}