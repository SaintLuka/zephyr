#include <stdexcept>
#include <cmath>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/matter/plast/plasticity.h>


namespace zephyr::phys {

Plasticity::Plasticity()
    : G(0.0), E(0.0), nu(0.5), Y(0.0) { }

Plasticity::Plasticity(const std::string &name) {
    table_params(name);
}

double Plasticity::shear() const {
    return G;
}

double Plasticity::young() const {
    return E;
}

double Plasticity::poisson() const {
    return nu;
}

double Plasticity::yield() const {
    return Y;
}

void Plasticity::table_params(const std::string& name) {
    if (name == "Al") {
        G  = 26.0_GPa;
        E  = 70.0_GPa;
        Y  = 18.0_MPa;
    }
    else if (name == "Fe") {
        G  = 52.5_GPa;
        E  = 180.0_GPa;
        Y  = 90.0_MPa;
    }
    else if (name == "Cu") {
        G  = 40.7_GPa;
        E  = 110.0_GPa;
        Y  = 33.0_MPa;
    }
    else {
        throw std::runtime_error("Unknown plasticity for '" + std::string(name) + "'");
    }

    nu = E / (2.0 * G) - 1.0;
}

} // namespace zephyr::phys
