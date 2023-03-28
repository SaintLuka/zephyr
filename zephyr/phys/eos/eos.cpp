#include <stdexcept>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>

#include <zephyr/phys/eos/eos.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/eos/stiffened_gas.h>

namespace zephyr { namespace phys {

double Eos::density_pe(double pressure, double energy) const {
    throw std::runtime_error("density_pe is not implemented");
}

double Eos::density_pt(double pressure, double temperature) const {
    throw std::runtime_error("density_pt is not implemented");
}

dRdE Eos::pressure_re(double density, double energy) const {
    throw std::runtime_error("pressure_re is not implemented");
}

double Eos::pressure_rt(double density, double temperature) const {
    throw std::runtime_error("pressure_rt is not implemented");
}

double Eos::energy_rp(double density, double pressure) const {
    throw std::runtime_error("energy_rp is not implemented");
}

double Eos::energy_rt(double density, double temperature) const {
    throw std::runtime_error("energy_rt is not implemented");
}

double Eos::energy_pt(double pressure, double temperature) const {
    throw std::runtime_error("energy_pt is not implemented");
}

double Eos::temperature_rp(double density, double pressure) const {
    throw std::runtime_error("temperature_rp not implemented");
}

double Eos::temperature_re(double density, double energy) const {
    throw std::runtime_error("temperature_re is not implemented");
}

double Eos::temperature_pe(double pressure, double energy) const {
    throw std::runtime_error("temperature_pe is not implemented");
}

double Eos::sound_speed_rp(double density, double pressure) const {
    throw std::runtime_error("sound_speed_rp is not implemented");
}

double Eos::sound_speed_re(double density, double energy) const {
    throw std::runtime_error("sound_speed_re is not implemented");
}

double Eos::sound_speed_rt(double density, double temperature) const {
    throw std::runtime_error("sound_speed_rt is not implemented");
}

double Eos::sound_speed_pe(double pressure, double energy) const {
    throw std::runtime_error("sound_speed_pe is not implemented");
}

double Eos::sound_speed_pt(double pressure, double temperature) const {
    throw std::runtime_error("sound_speed_pt is not implemented");
}

void test_eos(Eos& eos, double temperature, double density) {
    double pressure = eos.pressure_rt(density, temperature);
    double energy = eos.energy_rt(density, temperature);

    std::vector<double> err1 = {
            std::abs(density - eos.density_pe(pressure, energy)) / std::abs(density),
            std::abs(density - eos.density_pt(pressure, temperature)) / std::abs(density),
            std::abs(pressure - eos.pressure_re(density, energy)) / std::abs(pressure),
            std::abs(pressure - eos.pressure_rt(density, temperature)) / std::abs(pressure),
            std::abs(energy - eos.energy_rp(density, pressure)) / std::abs(energy),
            std::abs(energy - eos.energy_rt(density, temperature)) / std::abs(energy),
            std::abs(energy - eos.energy_pt(pressure, temperature)) / std::abs(energy),
            std::abs(temperature - eos.temperature_rp(density, pressure)) / std::abs(temperature),
            std::abs(temperature - eos.temperature_re(density, energy)) / std::abs(temperature),
            std::abs(temperature - eos.temperature_pe(pressure, energy)) / std::abs(temperature)
    };

    if (*std::max_element(err1.begin(), err1.end()) > 1.0e-14) {
        throw std::runtime_error("Eos test failed #1");
    }

    auto P = eos.pressure_re(density, energy);
    double c = std::sqrt(P.dR + P.value * P.dE / density / density);

    std::vector<double> err2 = {
            std::abs(c - eos.sound_speed_rp(density, pressure)) / c,
            std::abs(c - eos.sound_speed_re(density, energy)) / c,
            std::abs(c - eos.sound_speed_rt(density, temperature)) / c,
            std::abs(c - eos.sound_speed_pe(pressure, energy)) / c,
            std::abs(c - eos.sound_speed_pt(pressure, temperature)) / c
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
