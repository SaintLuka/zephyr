#include <stdexcept>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

#include <zephyr/phys/eos/eos.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/eos/stiffened_gas.h>

namespace zephyr { namespace phys {

dRdE Eos::pressure_re(double rho, double e, const Options& options) const {
    throw std::runtime_error("pressure_re is not implemented");
}

dRdP Eos::energy_rp(double rho, double p, const Options& options) const {
    throw std::runtime_error("energy_rp is not implemented");
}

double Eos::sound_speed_re(double rho, double e, const Options& options) const {
    auto opts = options;
    opts.deriv = true;

    auto P = pressure_re(rho, e, opts);
    double c2 = P.dR + P.value * P.dE / (rho * rho);
    return std::sqrt(c2);
}

double Eos::sound_speed_rp(double rho, double p, const Options& options) const {
    auto opts = options;
    opts.deriv = true;

    auto E = energy_rp(rho, p, opts);
    double c2 = (p / (rho * rho) - E.dR) / E.dP;
    return std::sqrt(c2);
}

double Eos::pressure_rt(double rho, double T, const Options& options) const {
    throw std::runtime_error("pressure_rt is not implemented");
}

dPdT Eos::density_pt(double p, double T, const Options& options) const {
    throw std::runtime_error("density_pt is not implemented");
}

dPdT Eos::energy_pt(double p, double T, const Options& options) const {
    throw std::runtime_error("energy_pt is not implemented");
}

double Eos::temperature_rp(double rho, double p, const Options& options) const {
    throw std::runtime_error("temperature_rp is not implemented");
}

double Eos::stiff_gamma(double rho, double p, const Options& options) const {
    throw std::runtime_error("stiff_gamma is not implemented");
}

double Eos::stiff_p0(double rho, double p, const Options& options) const {
    throw std::runtime_error("stiff_p0 is not implemented");
}

double Eos::stiff_e0(double rho, double p, const Options& options) const {
    throw std::runtime_error("stiff_e0 is not implemented");
}

void test_eos(Eos& eos, double T, double rho) {
    double p = eos.pressure_rt(rho, T);
    double e = eos.energy_rp(rho, p);

    std::vector<double> err1 = {
            std::abs(rho - eos.density_pt(p, T)) / std::abs(rho),
            std::abs(p - eos.pressure_re(rho, e)) / std::abs(p),
            std::abs(p - eos.pressure_rt(rho, T)) / std::abs(p),
            std::abs(e - eos.energy_rp(rho, p)) / std::abs(e),
            std::abs(e - eos.energy_pt(p, T)) / std::abs(e),
            std::abs(T - eos.temperature_rp(rho, p)) / std::abs(T),
    };

    if (*std::max_element(err1.begin(), err1.end()) > 1.0e-14) {
        throw std::runtime_error("Eos test failed #1");
    }

    auto P = eos.pressure_re(rho, e, {.deriv = true});
    double c = std::sqrt(P.dR + P.value * P.dE / (rho * rho));

    std::vector<double> err2 = {
            std::abs(c - eos.sound_speed_rp(rho, p)) / c,
            std::abs(c - eos.sound_speed_re(rho, e)) / c
    };

    if (*std::max_element(err2.begin(), err2.end()) > 1.0e-14) {
        throw std::runtime_error("Eos test failed #2");
    }
}

void Eos::test() {
    IdealGas gas1("Air");
    test_eos(gas1, 10.0_C, 1.3_kg_m3);

    StiffenedGas gas2("Water2");
    test_eos(gas2, 10.0_C, 1.1_g_cm3);
}

}
}
