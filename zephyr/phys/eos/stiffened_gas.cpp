#include <string>
#include <stdexcept>
#include <cmath>

#include <zephyr/phys/eos/stiffened_gas.h>

namespace zephyr { namespace phys {

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

double StiffenedGas::density_pe(double pressure, double energy) const {
    return (pressure + gamma * p_inf) / ((gamma - 1.0) * (energy - eps_0));
}

double StiffenedGas::density_pt(double pressure, double temperature) const {
    return (pressure + p_inf) / ((gamma - 1.0) * Cv * temperature);
}

dRdE StiffenedGas::pressure_re(double density, double energy) const {
    return {
            (gamma - 1.0) * density * (energy - eps_0) - gamma * p_inf,
            (gamma - 1.0) * (energy - eps_0),
            (gamma - 1.0) * density
    };
}

double StiffenedGas::pressure_rt(double density, double temperature) const {
    return (gamma - 1.0) * Cv * density * temperature  - p_inf;
}

double StiffenedGas::energy_rp(double density, double pressure) const {
    return eps_0 + (pressure + gamma * p_inf) / ((gamma - 1.0) * density);
}

double StiffenedGas::energy_rt(double density, double temperature) const {
    return Cv * temperature + eps_0 + p_inf / density;
}

double StiffenedGas::energy_pt(double pressure, double temperature) const {
    return eps_0 + ((pressure + gamma * p_inf) / (pressure + p_inf)) * Cv * temperature;
}

double StiffenedGas::temperature_rp(double density, double pressure) const {
    return (pressure + p_inf) / ((gamma - 1.0) * Cv * density);
}

double StiffenedGas::temperature_re(double density, double energy) const {
    return ((density * (energy - eps_0) - p_inf) / density) / Cv;
}

double StiffenedGas::temperature_pe(double pressure, double energy) const {
    return ((energy - eps_0) * (pressure + p_inf)) / (Cv * (pressure + gamma * p_inf));
}

double StiffenedGas::sound_speed_rp(double density, double pressure) const {
    return std::sqrt(gamma * (pressure + p_inf) / density);
}

double StiffenedGas::sound_speed_re(double density, double energy) const {
    return std::sqrt(gamma * (gamma - 1.0) * ((energy - eps_0) - p_inf / density));
}

double StiffenedGas::sound_speed_rt(double density, double temperature) const {
    return std::sqrt(gamma * (gamma - 1.0) * Cv * temperature);
}

double StiffenedGas::sound_speed_pe(double pressure, double energy) const {
    return std::sqrt(gamma * (gamma - 1.0) * (energy - eps_0) *
                     (pressure + p_inf) / (pressure + gamma * p_inf));
}

double StiffenedGas::sound_speed_pt(double pressure, double temperature) const {
    return std::sqrt(gamma * (gamma - 1.0) * Cv * temperature);
}

}
}