#include <string>
#include <stdexcept>
#include <cmath>

#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr { namespace phys {

TaitEos::TaitEos(double gamma, double Cv)
    : gamma(gamma), Cv(Cv) { }

TaitEos::TaitEos(const std::string &name) {
    if (name == "Air") {
        gamma = 1.4;
        Cv = 718.0_J_kgK;
    } else {
        throw std::runtime_error("Unknown Tait material '" + name + "'");
    }
}

dRdE TaitEos::pressure_re(double rho, double e, const Options& options) const {
    dRdE res {(gamma - 1.0) * rho * e };
    if (options.deriv) {
        res.dR = (gamma - 1.0) * e;
        res.dE = (gamma - 1.0) * rho;
    }
    return res;
}

dRdP TaitEos::energy_rp(double rho, double p, const Options& options) const {
    dRdP res { p / ((gamma - 1.0) * rho) };
    if (options.deriv) {
        res.dR = - res.value / rho;
        res.dP = 1.0 / ((gamma - 1.0) * rho);
    }
    return res;
}

double TaitEos::sound_speed_re(double rho, double e, const Options& options) const {
    return std::sqrt(gamma * (gamma - 1.0) * e);
}

double TaitEos::sound_speed_rp(double rho, double p, const Options& options) const {
    return std::sqrt(gamma * p / rho);
}

double TaitEos::pressure_rt(double rho, double T, const Options& options) const {
    return (gamma - 1.0) * Cv * rho * T;
}

dPdT TaitEos::density_pt(double p, double T, const Options& options) const {
    dPdT res{p / ((gamma - 1.0) * Cv * T)};
    if (options.deriv) {
        res.dP = 1.0 / ((gamma - 1.0) * Cv * T);
        res.dT = -res.value / T;
    }
    return res;
}

inline double sqr(double x) { return x * x; }

dPdT TaitEos::energy_pt(double p, double T, const Options& options) const {
    dPdT res{Cv * T};
    if (options.deriv) {
        res.dP = 0.0;
        res.dT = Cv;
    }
    return res;
}

double TaitEos::temperature_rp(double rho, double p, const Options& options) const {
    return p / ((gamma - 1.0) * Cv * rho);
}

double TaitEos::stiff_gamma(double rho, double p, const Options& options) const {
    return gamma;
}

double TaitEos::stiff_p0(double rho, double p, const Options& options) const {
    return 0.0;
}

double TaitEos::stiff_e0(double rho, double p, const Options& options) const {
    return 0.0;
}

}
}