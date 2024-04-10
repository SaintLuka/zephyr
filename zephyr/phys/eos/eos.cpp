#include <stdexcept>
#include <cmath>

#include <zephyr/phys/eos/eos.h>
#include <zephyr/phys/eos/stiffened_gas.h>

namespace zephyr { namespace phys {

dRdE Eos::pressure_re(double rho, double eps, const Options& options) const {
    throw std::runtime_error("pressure_re is not implemented");
}

double Eos::energy_rp(double rho, double P, const Options& options) const {
    throw std::runtime_error("energy_rp is not implemented");
}

double Eos::sound_speed_re(double rho, double eps, const Options& options) const {
    auto opts = options;
    opts.deriv = true;

    auto P = pressure_re(rho, eps, opts);
    double c2 = P.dR + P.val * P.dE / (rho * rho);
    return std::sqrt(c2);
}

double Eos::sound_speed_rp(double rho, double P, const Options& options) const {
    double eps = energy_rp(rho, P, options);
    return sound_speed_re(rho, eps);
}

double Eos::pressure_rt(double rho, double T, const Options& options) const {
    throw std::runtime_error("pressure_rt is not implemented");
}

double Eos::temperature_rp(double rho, double P, const Options& options) const {
    throw std::runtime_error("temperature_rp is not implemented");
}

dPdT Eos::volume_pt(double P, double T, const Options& options) const {
    throw std::runtime_error("volume_pt is not implemented");
}

dPdT Eos::energy_pt(double P, double T, const Options& options) const {
    throw std::runtime_error("energy_pt is not implemented");
}

StiffenedGas Eos::stiffened_gas(double rho, double P, const Options& options) const {
    double eps = energy_rp(rho, P, options);
    dRdE Pre = pressure_re(rho, eps, {.deriv = true});

    double gamma = 1.0 + Pre.dE / rho;
    double eps_0 = eps - rho * Pre.dR / Pre.dE;
    double P0 = (rho * Pre.dR - P) / gamma;

    return StiffenedGas(gamma, P0, eps_0, NAN);
}

double Eos::min_pressure() const {
    throw std::runtime_error("min_pressure is not implemented");
}

}
}
