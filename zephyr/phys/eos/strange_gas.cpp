#include <string>
#include <stdexcept>
#include <cmath>

#include <zephyr/phys/eos/strange_gas.h>

namespace zephyr::phys {

StrangeGas::StrangeGas(double gamma, double p_inf, double eps_0, double Cv, double T0)
    : gamma(gamma), P0(p_inf), eps_0(eps_0), Cv(Cv), T0(T0) { }

StrangeGas::StrangeGas(const std::string &name) {
    // @formatter:off
    T0 = 0.0;
    Pe = 100000.0;

    if (name == "Air") {
        gamma = 1.4;
        Cv    = 718.0_J_kgK;
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
        T0    = 221.3;
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

double StrangeGas::P_fix(double P) const {
    return 0.5 * (P + std::sqrt(P*P + Pe*Pe));
}
    
double StrangeGas::dP_fix(double P) const {
    return 0.5 * (P / std::sqrt(P * P + Pe * Pe) + 1);
}

double StrangeGas::P_inv(double P) const {
    return (4.0 * P * P - Pe * Pe) / (4.0 * P);
}

double StrangeGas::dP_inv(double P) const {
    return Pe * Pe / (4.0 * P * P) + 1.0;
}


dRdE StrangeGas::pressure_re(double rho, double eps, const Options& options) const {
    double Pl = (gamma - 1.0) * rho * (eps - eps_0) - gamma * P0;
    dRdE res {P_inv(Pl)};
    if (options.deriv) {
        res.dR = dP_inv(Pl) * (gamma - 1.0) * (eps - eps_0);
        res.dE = dP_inv(Pl) * (gamma - 1.0) * rho;
    }
    return res;
}

double StrangeGas::energy_rp(double rho, double P, const Options& options) const {
    return eps_0 + (P_fix(P) + gamma * P0) / ((gamma - 1.0) * rho);
}
/*
double StrangeGas::sound_speed_re(double rho, double eps, const Options& options) const {
    return std::sqrt(gamma * (gamma - 1.0) * ((eps - eps_0) - P0 / rho));
}

double StrangeGas::sound_speed_rp(double rho, double P, const Options& options) const {
    return std::sqrt(gamma * (P + P0) / rho);
}
*/
double StrangeGas::pressure_rt(double rho, double T, const Options& options) const {
    return P_inv((gamma - 1.0) * Cv * rho * (T - T0) - P0);
}

double StrangeGas::temperature_rp(double rho, double P, const Options& options) const {
    return T0 + (P_fix(P) + P0) / ((gamma - 1.0) * Cv * rho);
}

dPdT StrangeGas::volume_pt(double P, double T, const Options& options) const {
    dPdT res{((gamma - 1.0) * Cv * (T - T0)) / (P_fix(P) + P0)};
    if (options.deriv) {
        res.dP = -res.val / (P_fix(P) + P0) * dP_fix(P);
        res.dT = +res.val / (T - T0);
    }
    return res;
}

inline double sqr(double x) { return x * x; }

dPdT StrangeGas::energy_pt(double P, double T, const Options& options) const {
    dPdT res{eps_0 + ((P_fix(P) + gamma * P0) / (P_fix(P) + P0)) * Cv * (T - T0)};
    if (options.deriv) {
        res.dP = (1.0 - gamma) * P0 * Cv * (T - T0) / sqr(P_fix(P) + P0) * dP_fix(P);
        res.dT = ((P_fix(P) + gamma * P0) / (P_fix(P) + P0)) * Cv;
    }
    return res;
}

/*
StiffenedGas StrangeGas::stiffened_gas(double rho, double P, const Options& options) const {
    return StiffenedGas(gamma, P0, eps_0, Cv, T0);
}
*/
double StrangeGas::min_pressure() const {
    return -1.0e10;
}

}
