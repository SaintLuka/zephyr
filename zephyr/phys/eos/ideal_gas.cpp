#include <string>
#include <stdexcept>
#include <cmath>

#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/eos/stiffened_gas.h>

namespace zephyr { namespace phys {

IdealGas::IdealGas(double gamma, double Cv)
    : gamma(gamma), Cv(Cv) { }

IdealGas::IdealGas(const std::string &name) {
    if (name == "Air") {
        gamma = 1.4;
        Cv = 718.0_J_kgK;
    } else {
        throw std::runtime_error("Unknown Ideal gas '" + name + "'");
    }
}

dRdE IdealGas::pressure_re(double rho, double eps, const Options& options) const {
    dRdE res {(gamma - 1.0) * rho * eps };
    if (options.deriv) {
        res.dR = (gamma - 1.0) * eps;
        res.dE = (gamma - 1.0) * rho;
    }
    return res;
}

double IdealGas::energy_rp(double rho, double P, const Options& options) const {
    return P / ((gamma - 1.0) * rho);
}

double IdealGas::sound_speed_re(double rho, double eps, const Options& options) const {
    return std::sqrt(gamma * (gamma - 1.0) * eps);
}

double IdealGas::sound_speed_rp(double rho, double P, const Options& options) const {
    return std::sqrt(gamma * P / rho);
}

double IdealGas::pressure_rt(double rho, double T, const Options& options) const {
    return (gamma - 1.0) * Cv * rho * T;
}

double IdealGas::temperature_rp(double rho, double P, const Options& options) const {
    return P / ((gamma - 1.0) * Cv * rho);
}

dPdT IdealGas::volume_pt(double P, double T, const Options& options) const {
    dPdT res{((gamma - 1.0) * Cv * T) / P};
    if (options.deriv) {
        res.dP = -res.val / P;
        res.dT = +res.val / T;
    }
    return res;
}

dPdT IdealGas::energy_pt(double P, double T, const Options& options) const {
    dPdT res{Cv * T};
    if (options.deriv) {
        res.dP = 0.0;
        res.dT = Cv;
    }
    return res;
}

StiffenedGas IdealGas::stiffened_gas(double rho, double p, const Options& options) const {
    return StiffenedGas(gamma, 0.0, 0.0, Cv);
}

}
}