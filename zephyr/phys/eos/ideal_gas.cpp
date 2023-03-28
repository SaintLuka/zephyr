#include <string>
#include <stdexcept>
#include <cmath>

#include <zephyr/phys/eos/ideal_gas.h>

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

double IdealGas::density_pe(double pressure, double energy) const {
    return pressure / ((gamma - 1.0) * energy);
}

double IdealGas::density_pt(double pressure, double temperature) const {
    return pressure / ((gamma - 1.0) * Cv * temperature);
}

dRdE IdealGas::pressure_re(double density, double energy) const {
    return {
            (gamma - 1.0) * density * energy,
            (gamma - 1.0) * energy,
            (gamma - 1.0) * density
    };
}

double IdealGas::pressure_rt(double density, double temperature) const {
    return (gamma - 1.0) * Cv * density * temperature;
}

double IdealGas::energy_rp(double density, double pressure) const {
    return pressure / ((gamma - 1.0) * density);
}

double IdealGas::energy_rt(double density, double temperature) const {
    return Cv * temperature;
}

double IdealGas::energy_pt(double pressure, double temperature) const {
    return Cv * temperature;
}

double IdealGas::temperature_rp(double density, double pressure) const {
    return pressure / ((gamma - 1.0) * Cv * density);
}

double IdealGas::temperature_re(double density, double energy) const {
    return energy / Cv;
}

double IdealGas::temperature_pe(double pressure, double energy) const {
    return energy / Cv;
}

double IdealGas::sound_speed_rp(double density, double pressure) const {
    return std::sqrt(gamma * pressure / density);
}

double IdealGas::sound_speed_re(double density, double energy) const {
    return std::sqrt(gamma * (gamma - 1.0) * energy);
}

double IdealGas::sound_speed_rt(double density, double temperature) const {
    return std::sqrt(gamma * (gamma - 1.0) * Cv * temperature);
}

double IdealGas::sound_speed_pe(double pressure, double energy) const {
    return std::sqrt(gamma * (gamma - 1.0) * energy);
}

double IdealGas::sound_speed_pt(double pressure, double temperature) const {
    return std::sqrt(gamma * (gamma - 1.0) * Cv * temperature);
}

}
}