#include <zephyr/phys/tests/ivp.h>

namespace zephyr::phys {

using mesh::EuCell;

const Materials& IVP::materials() const {
    return m_materials;
}

MixturePT IVP::mixture_PT() const {
    return m_materials.mixture_PT();
}

Eos::Ptr IVP::get_eos(int idx) const {
    assert(idx < m_materials.size());
    return m_materials[idx].eos();
}

double IVP::energy(const Vector3d &r) const {
    return m_materials[index(r)].energy_rP(density(r), pressure(r));
}

double IVP::temperature(const Vector3d &r) const {
    return m_materials[index(r)].temperature_rP(density(r), pressure(r));
}

Fractions IVP::fractions(const Vector3d &r) const {
    return Fractions::Pure(index(r));
}

bool IVP::inside(const Vector3d &r, int idx) const {
    return idx == 0;
}

int IVP::index_t(const Vector3d& r, double t) const {
    return 0;
}

double IVP::density_t(const Vector3d &r, double t) const {
    return NAN;
}

Vector3d IVP::velocity_t(const Vector3d &r, double t) const {
    return {NAN, NAN, NAN};
}

double IVP::pressure_t(const Vector3d &r, double t) const {
    return NAN;
}

double IVP::energy_t(const Vector3d &r, double t) const {
    double rho = density_t(r, t);
    double P = pressure_t(r, t);
    int idx = index_t(r, t);
    if (!std::isnan(rho) && !std::isnan(P) && idx < m_materials.size()) {
        return m_materials[index_t(r, t)].energy_rP(rho, P);
    }
    return NAN;
}

double IVP::temperature_t(const Vector3d& r, double t) const {
    double rho = density_t(r, t);
    double P = pressure_t(r, t);
    int idx = index_t(r, t);
    if (!std::isnan(rho) && !std::isnan(P) && idx < m_materials.size()) {
        return m_materials[index_t(r, t)].temperature_rP(rho, P);
    }
    return NAN;
}

Fractions IVP::fractions_t(const Vector3d &r, double t) const {
    return Fractions::Pure(index_t(r, t));
}

bool IVP::inside_t(const Vector3d &r, int idx, double t) const {
    return idx == 0;
}

double IVP::density_mean(EuCell& cell, int n) const {
    auto get_density = [this](const Vector3d &r) -> double {
        return this->density(r);
    };

    if (n < 2 || cell.const_function(get_density)) {
        return density(cell.center());
    }

    double V = cell.volume();
    double M = cell.integrate_low(get_density, n);
    return M / V;
}

Vector3d IVP::momentum_mean(EuCell& cell, int n) const {
    auto momentum_x = [this](const Vector3d &r) -> double {
        return this->density(r) * this->velocity(r).x();
    };
    auto momentum_y = [this](const Vector3d &r) -> double {
        return this->density(r) * this->velocity(r).y();
    };
    auto momentum_z = [this](const Vector3d &r) -> double {
        return this->density(r) * this->velocity(r).z();
    };

    if (n < 2 ||
        cell.const_function(momentum_x) ||
        cell.const_function(momentum_y) ||
        cell.const_function(momentum_z)) {

        return {momentum_x(cell.center()),
                momentum_y(cell.center()),
                momentum_z(cell.center())};
        }

    double V = cell.volume();
    Vector3d P = {
            cell.integrate_low(momentum_x, n),
            cell.integrate_low(momentum_y, n),
            cell.integrate_low(momentum_z, n)
    };
    return P / V;
}

double IVP::energy_mean(EuCell& cell, int n) const {
    auto get_energy = [this](const Vector3d &r) -> double {
        Vector3d v = this->velocity(r);
        return this->density(r) * (this->energy(r) + 0.5 * v.dot(v));
    };

    if (n < 2 || cell.const_function(get_energy)) {
        return get_energy(cell.center());
    }

    double V = cell.volume();
    double E = cell.integrate_low(get_energy, n);
    return E / V;
}

Fractions IVP::mass_fractions(EuCell& cell, int n) const {
    throw std::runtime_error("IVP::masses");
}

Fractions IVP::volume_fractions(EuCell& cell, int n) const {
    throw std::runtime_error("IVP::volumes");
}

} // namespace zephyr::phys