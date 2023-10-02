#include <string>
#include <stdexcept>
#include <cmath>

#include <zephyr/phys/eos/stiffened_gas.h>

namespace zephyr { namespace phys {

StiffenedGas::StiffenedGas(double gamma, double p_inf, double eps_0, double Cv)
    : gamma(gamma), p_inf(p_inf), eps_0(eps_0), Cv(Cv) { }

StiffenedGas::StiffenedGas(const std::string &name) {
    // @formatter:off

    if (name == "Air") {
        gamma = 1.4;
        Cv    = 0.0;
        p_inf = 0.0;
        eps_0 = 0.0;
    }
    else if (name == "SF6") {
        gamma = 1.076;
        Cv    = 96.6_J_kgK;
        p_inf = 0.0;
        eps_0 = 0.0;
    }
    else if (name == "Water") {
        gamma = 4.4;
        Cv    = 1400.0_J_kgK;
        p_inf = 600.0_MPa;
        eps_0 = 0.0;
    }
    else if (name == "Water2") {
        gamma = 3.0;
        Cv    = 1400.0_J_kgK;
        p_inf = 853.3_MPa;
        eps_0 = -1.148_MJ_kg;
    }
    else if (name == "Copper") {
        gamma = 4.0;
        Cv    = 384.0_J_kgK;
        p_inf = 341.0;
        eps_0 = 0.0;
    }
    else {
        throw std::runtime_error("Unknown stiffened gas '" + std::string(name) + "'");
    }
    // @formatter:on
}

dRdE StiffenedGas::pressure_re(double rho, double e, const Options& options) const {
    dRdE res {(gamma - 1.0) * rho * (e - eps_0) - gamma * p_inf };
    if (options.deriv) {
        res.dR = (gamma - 1.0) * (e - eps_0);
        res.dE = (gamma - 1.0) * rho;
    }
    return res;
}

dRdP StiffenedGas::energy_rp(double rho, double p, const Options& options) const {
    dRdP res { eps_0 + (p + gamma * p_inf) / ((gamma - 1.0) * rho) };
    if (options.deriv) {
        res.dR = (eps_0 - res.value) / rho;
        res.dP = 1.0 / ((gamma - 1.0) * rho);
    }
    return res;
}

double StiffenedGas::sound_speed_re(double rho, double e, const Options& options) const {
    return std::sqrt(gamma * (gamma - 1.0) * ((e - eps_0) - p_inf / rho));
}

double StiffenedGas::sound_speed_rp(double rho, double p, const Options& options) const {
    return std::sqrt(gamma * (p + p_inf) / rho);
}

double StiffenedGas::pressure_rt(double rho, double T, const Options& options) const {
    return (gamma - 1.0) * Cv * rho * T - p_inf;
}

dPdT StiffenedGas::density_pt(double p, double T, const Options& options) const {
    dPdT res{(p + p_inf) / ((gamma - 1.0) * Cv * T)};
    if (options.deriv) {
        res.dP = 1.0 / ((gamma - 1.0) * Cv * T);
        res.dT = -res.value / T;
    }
    return res;
}

inline double sqr(double x) { return x * x; }

dPdT StiffenedGas::energy_pt(double p, double T, const Options& options) const {
    dPdT res{eps_0 + ((p + gamma * p_inf) / (p + p_inf)) * Cv * T};
    if (options.deriv) {
        res.dP = (1.0 - gamma) * p_inf * Cv * T / sqr(p - p_inf);
        res.dT = ((p + gamma * p_inf) / (p + p_inf)) * Cv;
    }
    return res;
}

double StiffenedGas::temperature_rp(double rho, double p, const Options& options) const {
    return (p + p_inf) / ((gamma - 1.0) * Cv * rho);
}

double StiffenedGas::stiff_gamma(double rho, double p, const Options& options) const {
    return gamma;
}

double StiffenedGas::stiff_p0(double rho, double p, const Options& options) const {
    return p_inf;
}

double StiffenedGas::stiff_e0(double rho, double p, const Options& options) const {
    return eps_0;
}

}
}