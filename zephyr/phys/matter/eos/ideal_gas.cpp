#include <string>
#include <stdexcept>
#include <cmath>

#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/matter/eos/stiffened_gas.h>

namespace zephyr::phys {

IdealGas::IdealGas(double gamma, double Cv)
    : gamma(gamma), Cv(Cv) {

    rho_0 = 1.0 / volume_PT(1.0e5, 300);
}

IdealGas::IdealGas(const std::string &name) {
    if (name == "Air") {
        gamma = 1.4;
        Cv    = 718.0_J_kgK;
        rho_0 = 1.25_kg_m3;
    }
    else if (name == "SF6") {
        gamma = 1.1074;
        Cv    = 567.0_J_kgK;
        rho_0 = 6.17_kg_m3;
    } else {
        throw std::runtime_error("Unknown Ideal gas '" + name + "'");
    }
}

dRdE IdealGas::pressure_re(double rho, double e, const Options& options) const {
    dRdE res {(gamma - 1.0) * rho * e };
    if (options.deriv) {
        res.dR = (gamma - 1.0) * e;
        res.dE = (gamma - 1.0) * rho;
    }
    return res;
}

dRdT IdealGas::pressure_rT(double rho, double T, const Options& options) const {
    dRdT P = (gamma - 1.0) * Cv * rho * T;
    if (options.deriv) {
        P.dR = (gamma - 1.0) * Cv * T;
        P.dT = (gamma - 1.0) * Cv * rho;
    }
    return P;
}

dRdT IdealGas::energy_rT(double rho, double T, const Options& options) const {
    dRdT e = Cv * T;
    if (options.deriv) {
        e.dR = 0.0;
        e.dT = Cv;
    }
    return e;
}

double IdealGas::sound_speed_re(double rho, double e, const Options& options) const {
    return std::sqrt(gamma * (gamma - 1.0) * e);
}

double IdealGas::sound_speed_rP(double rho, double P, const Options& options) const {
    return std::sqrt(gamma * P / rho);
}

double IdealGas::energy_rP(double rho, double P, const Options& options) const {
    return P / ((gamma - 1.0) * rho);
}

double IdealGas::temperature_rP(double rho, double P, const Options& options) const {
    return P / ((gamma - 1.0) * Cv * rho);
}

dPdT IdealGas::volume_PT(double P, double T, const Options& options) const {
    dPdT res{((gamma - 1.0) * Cv * T) / P};
    if (options.deriv) {
        res.dP = -res.val / P;
        res.dT = +res.val / T;
    }
    return res;
}

dPdT IdealGas::energy_PT(double P, double T, const Options& options) const {
    dPdT res{Cv * T};
    if (options.deriv) {
        res.dP = 0.0;
        res.dT = Cv;
    }
    return res;
}

StiffenedGas IdealGas::stiffened_gas(double rho, double P, const Options& options) const {
    return StiffenedGas(gamma, 0.0, 0.0, Cv);
}

double IdealGas::min_pressure() const {
    return 0.0;
}

void IdealGas::adjust_cv(double rho, double P, double T) {
    Cv = P / ((gamma - 1.0) * Cv * rho);
}

} // namespace zephyr::phys