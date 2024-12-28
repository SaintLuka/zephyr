#include <stdexcept>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/matter/cond/conductivity.h>


namespace zephyr::phys {

Conductivity::Conductivity(double kappa)
    : m_kappa(kappa) {
}

Conductivity::Conductivity(const std::string& name) {
    table_params(name);
}

double Conductivity::kappa(double temperature) const {
    return m_kappa;
}

void Conductivity::table_params(const std::string& name) {
    if (name == "Air") {
        m_kappa = 0.022;
    }
    else if (name == "Water") {
        m_kappa = 0.6_W_mK;
    }
    else if (name == "Cu") {
        m_kappa = 401.0_W_mK;
    }
    else if (name == "Al") {
        m_kappa = 220.0_W_mK;
    }
    else if (name == "Fe") {
        m_kappa = 92.0_W_mK;
    }
    else if (name == "Pt") {
        m_kappa = 70.0_W_mK;
    }
    else {
        throw std::runtime_error("Unknown conductivity for '" + std::string(name) + "'");
    }
}

} // namespace zephyr::phys
