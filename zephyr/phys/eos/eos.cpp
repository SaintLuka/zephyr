#include <stdexcept>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

#include <zephyr/phys/eos/eos.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/eos/stiffened_gas.h>
#include <iostream>

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
    double c2 = P.dR + P.value * P.dE / (rho * rho);
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
    throw std::runtime_error("stiff_gamma is not implemented");
}

inline double error(double x, double y) {
    return std::abs((x - y) / (std::abs(x) + std::abs(y)));
}

void test_eos(Eos& eos, double T, double rho) {
    double p = eos.pressure_rt(rho, T);
    double e = eos.energy_rp(rho, p);

    std::vector<double> err1 = {
            error(p, eos.pressure_re(rho, e)),
            error(e, eos.energy_rp(rho, p)),
            error(p, eos.pressure_rt(rho, T)),
            error(rho, 1.0/eos.volume_pt(p, T)),
            error(e, eos.energy_pt(p, T)),
            error(T, eos.temperature_rp(rho, p))
    };

    if (*std::max_element(err1.begin(), err1.end()) > 1.0e-14) {
        throw std::runtime_error("Eos test failed #1");
    }

    auto P = eos.pressure_re(rho, e, {.deriv = true});
    auto V = eos.volume_pt(p, T, {.deriv = true});
    auto E = eos.energy_pt(p, T, {.deriv = true});

    double c = std::sqrt(P.dR + P.value * P.dE / (rho * rho));

    std::vector<double> err2 = {
            std::abs(c - eos.sound_speed_rp(rho, p)) / c,
            std::abs(c - eos.sound_speed_re(rho, e)) / c
    };

    if (*std::max_element(err2.begin(), err2.end()) > 1.0e-14) {
        throw std::runtime_error("Eos test failed #2");
    }

    // Проверка производных
    double eps = 1.0e-5;

    double dP = eps * p;
    double dR = eps * rho;
    double dT = eps * T;
    double dE = eps * e;

    double dPdR_e = (eos.pressure_re(rho + dR, e) - eos.pressure_re(rho - dR, e)) / (2 * dR);
    double dPdE_r = (eos.pressure_re(rho, e + dE) - eos.pressure_re(rho, e - dE)) / (2 * dE);
    double dVdP_t = (eos.volume_pt(p + dP, T) - eos.volume_pt(p - dP, T)) / (2 * dP);
    double dVdT_p = (eos.volume_pt(p, T + dT) - eos.volume_pt(p, T - dT)) / (2 * dT);
    double dEdP_t = (eos.energy_pt(p + dP, T) - eos.energy_pt(p - dP, T)) / (2*dP);
    double dEdT_p = (eos.energy_pt(p, T + dT) - eos.energy_pt(p, T - dT)) / (2*dT);

    std::vector<double> err3 = {
            error(P.dR, dPdR_e),
            error(P.dE, dPdE_r),
            error(V.dP, dVdP_t),
            error(V.dT, dVdT_p),
            error(E.dP, dEdP_t),
            error(E.dT, dEdT_p)
    };

    if (*std::max_element(err3.begin(), err3.end()) > 1.0e-8) {
        throw std::runtime_error("Eos test failed #3");
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
