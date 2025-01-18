#include <cmath>
#include <stdexcept>

#include <zephyr/phys/matter/eos/eos.h>
#include <zephyr/phys/matter/eos/stiffened_gas.h>

namespace zephyr::phys {

dRdE Eos::pressure_re(double rho, double e, const EosOptions& options) const {
    throw std::runtime_error("pressure_re is not implemented");
}

dRdT Eos::pressure_rT(double rho, double T, const EosOptions& options) const {
    throw std::runtime_error("pressure_rT is not implemented");
}

dRdP Eos::energy_rP(double rho, double P, const EosOptions& options) const {
    throw std::runtime_error("energy_rP is not implemented");
}

dRdT Eos::energy_rT(double rho, double T, const EosOptions& options) const {
    throw std::runtime_error("energy_rT is not implemented");
}

double Eos::temperature_rP(double rho, double P, const EosOptions& options) const {
    throw std::runtime_error("temperature_rP is not implemented");
}

double Eos::sound_speed_re(double rho, double e, const EosOptions& options) const {
    auto opts = options;
    opts.deriv = true;

    auto P = pressure_re(rho, e, opts);
    double c2 = P.dR + P.val * P.dE / (rho * rho);
    return std::sqrt(c2);
}

double Eos::sound_speed_rP(double rho, double P, const EosOptions& options) const {
    auto opts = options;
    opts.deriv = true;

    dRdP e = energy_rP(rho, P, opts);
    double c2 = (P / (rho * rho) - e.dR) / e.dP;
    return std::sqrt(c2);
}

dPdT Eos::volume_PT(double P, double T, const EosOptions& options) const {
    throw std::runtime_error("volume_PT is not implemented");
}

dPdT Eos::energy_PT(double P, double T, const EosOptions& options) const {
    throw std::runtime_error("energy_PT is not implemented");
}

StiffenedGas Eos::stiffened_gas(double rho, double P, const EosOptions& options) const {
    EosOptions opts=options;
    opts.deriv = true;
    dRdP e = energy_rP(rho, P, opts);
    double gamma = 1.0 + 1.0 / (rho * e.dP);
    double e0 = e + rho * e.dR;
    double P0 = -(rho * e.dR + P * e.dP) / (gamma * e.dP);
    return StiffenedGas(gamma, P0, e0, NAN);
}

double Eos::min_pressure() const {
    throw std::runtime_error("min_pressure is not implemented");
}

void Eos::adjust_cv(double rho, double P, double T) {
    throw std::runtime_error("adjust_cv is not implemented");
}

void Eos::adjust_T0(double rho, double P, double T) {
    throw std::runtime_error("adjust_T0 is not implemented");
}

} // namespace zephyr::phys
