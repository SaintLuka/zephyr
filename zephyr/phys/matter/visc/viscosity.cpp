#include <stdexcept>
#include <cmath>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/matter/visc/viscosity.h>


namespace zephyr::phys {

Viscosity::Viscosity()
    : nu(0.0), eta(0.0), zeta(0.0) { }

Viscosity::Viscosity(double nu)
    :nu(nu), eta(NAN), zeta(NAN) { }

Viscosity::Viscosity(double eta, double zeta)
    : nu(NAN), eta(eta), zeta(zeta) { }

Viscosity::Viscosity(const std::string &name) {
    table_params(name);
}

double Viscosity::kinematic_visc(double temperature) const {
    return nu;
}

double Viscosity::shear_visc(double temperature) const {
    return eta;
}

double Viscosity::volume_visc(double temperature) const {
    return zeta;
}

void Viscosity::table_params(const std::string& name) {
    if (name == "Water") {
        nu   = 8.9e-7;
        eta  = 8.9e-4;
        zeta = 0.0;
    }
    else {
        throw std::runtime_error("Unknown viscosity for '" + std::string(name) + "'");
    }
}

} // namespace zephyr::phys
