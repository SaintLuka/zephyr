#include <zephyr/phys/matter/material.h>
#include <zephyr/phys/matter/eos/stiffened_gas.h>

namespace zephyr::phys {

const std::string &Material::tag() const {
    return m_tag;
}

void Material::set_tag(const std::string &tag) {
    m_tag = tag;
}

void Material::add(Eos::Ref eos) {
    m_eos = eos;
}

void Material::add(Conductivity::Ref cond) {
    m_cond = cond;
}

void Material::add(Viscosity::Ref visc) {
    m_visc = visc;
}

void Material::add(Plasticity::Ref plast) {
    m_plast = plast;
}

void Material::operator+=(Eos::Ref eos) {
    m_eos = eos;
}

void Material::operator+=(Conductivity::Ref cond) {
    m_cond = cond;
}

void Material::operator+=(Viscosity::Ref visc) {
    m_visc = visc;
}

void Material::operator+=(Plasticity::Ref plast) {
    m_plast = plast;
}

double Material::density() const {
    return m_eos->density();
}

dRdE Material::pressure_re(double density, double energy, const EosOptions &options) const {
    return m_eos->pressure_re(density, energy, options);
}

dRdT Material::pressure_rT(double density, double temperature, const EosOptions &options) const {
    return m_eos->pressure_rT(density, temperature, options);
}

dRdT Material::energy_rT(double density, double temperature, const EosOptions &options) const {
    return m_eos->energy_rT(density, temperature, options);
}

double Material::sound_speed_re(double density, double energy, const EosOptions &options) const {
    return m_eos->sound_speed_re(density, energy, options);
}

double Material::sound_speed_rP(double density, double pressure, const EosOptions &options) const {
    return m_eos->sound_speed_rP(density, pressure, options);
}

double Material::energy_rP(double density, double pressure, const EosOptions &options) const {
    return m_eos->energy_rP(density, pressure, options);
}

double Material::temperature_rP(double density, double pressure, const EosOptions &options) const {
    return m_eos->temperature_rP(density, pressure, options);
}

dPdT Material::volume_PT(double pressure, double temperature, const EosOptions &options) const {
    return m_eos->volume_PT(pressure, temperature, options);
}

dPdT Material::energy_PT(double pressure, double temperature, const EosOptions &options) const {
    return m_eos->energy_PT(pressure, temperature, options);
}

StiffenedGas Material::stiffened_gas(double density, double pressure, const EosOptions &options) const {
    return m_eos->stiffened_gas(density, pressure, options);
}

double Material::min_pressure() const {
    return m_eos->min_pressure();
}

void Material::adjust_cv(double rho_ref, double P_ref, double T_ref) {
    return m_eos->adjust_cv(rho_ref, P_ref, T_ref);
}

void Material::adjust_T0(double rho_ref, double P_ref, double T_ref) {
    return m_eos->adjust_T0(rho_ref, P_ref, T_ref);
}

double Material::kappa(double temperature) const {
    return m_cond->kappa(temperature);
}

double Material::kinematic_visc(double temperature) const {
    return m_visc->kinematic_visc(temperature);
}

double Material::shear_visc(double temperature) const {
    return m_visc->shear_visc(temperature);
}

double Material::volume_visc(double temperature) const {
    return m_visc->volume_visc(temperature);
}

double Material::shear() const {
    return m_plast->shear();
}

double Material::young() const {
    return m_plast->young();
}

double Material::poisson() const {
    return m_plast->poisson();
}

double Material::yield() const {
    return m_plast->yield();
}

} // namespace zephyr::phys