#include <zephyr/math/cfd/models.h>
#include <boost/format.hpp>
#include <stdexcept>
#include <iostream>

namespace zephyr::math {

namespace smf {

PState::PState()
    : density(0.0),
      velocity({0.0, 0.0, 0.0}),
      pressure(0.0),
      energy(0.0)
{

}

PState::PState(const double &density, const Vector3d &velocity,
               const double &pressure, const double &energy)
        : density(density), velocity(velocity),
          pressure(pressure), energy(energy) {
}

PState::PState(const QState &q, const phys::Eos &eos) {
    density = q.mass;
    velocity = q.momentum / density;
    energy = q.energy / density - 0.5 * velocity.squaredNorm();
    pressure = eos.pressure_re(density, energy);
}

void PState::to_local(const Vector3d &normal) {
    Rotate::to_local(velocity, normal);
}

PState PState::in_local(const Vector3d &normal) const {
    PState z(*this);
    z.to_local(normal);
    return z;
}

void PState::to_global(const Vector3d &normal) {
    Rotate::to_global(velocity, normal);
}

PState PState::in_global(const Vector3d &normal) const {
    PState z(*this);
    z.to_global(normal);
    return z;
}

std::ostream &operator<<(std::ostream &os, const PState &state) {
    os << "density: " << state.density << " velocity: {" << state.velocity.x() << ", " << state.velocity.y() << ", "
       << state.velocity.z() <<
       "} pressure: " << state.pressure << " energy: " << state.energy;
    return os;
}

QState::QState(const double &mass, const Vector3d &momentum, const double &energy)
    : mass(mass),
      momentum(momentum),
      energy(energy) {
}

QState::QState(const PState &z) {
    mass = z.density;
    momentum = z.density * z.velocity;
    energy = z.density * (z.energy + 0.5 * z.velocity.squaredNorm());
}

void QState::to_local(const Vector3d &normal) {
    Rotate::to_local(momentum, normal);
}

QState QState::in_local(const Vector3d &normal) const {
    QState q(*this);
    q.to_local(normal);
    return q;
}

void QState::to_global(const Vector3d &normal) {
    Rotate::to_global(momentum, normal);
}

QState QState::in_global(const Vector3d &normal) const {
    QState q(*this);
    q.to_global(normal);
    return q;
}

std::ostream &operator<<(std::ostream &os, const QState &state) {
    os << "mass: " << state.mass <<
       " momentum: {" << state.momentum.x() << ", " << state.momentum.y() << ", " << state.momentum.z() <<
       "} energy: " << state.energy;
    return os;
}

Flux::Flux()
    : mass(0.0),
      momentum({0.0, 0.0, 0.0}),
      energy(0.0) {

}

Flux::Flux(double mass, const Vector3d &momentum, double energy)
    : mass(mass),
      momentum(momentum),
      energy(energy) {

}

Flux::Flux(const PState &z) {
    mass = z.density * z.velocity.x();
    momentum.x() = z.density * z.velocity.x() * z.velocity.x() + z.pressure;
    momentum.y() = z.density * z.velocity.x() * z.velocity.y();
    momentum.z() = z.density * z.velocity.x() * z.velocity.z();
    energy = (z.density * (z.energy + 0.5 * z.velocity.squaredNorm()) + z.pressure) * z.velocity.x();
}

void Flux::to_local(const Vector3d &normal) {
    Rotate::to_local(momentum, normal);
}

Flux Flux::in_local(const Vector3d &normal) const {
    Flux f(*this);
    f.to_local(normal);
    return f;
}

void Flux::to_global(const Vector3d &normal) {
    Rotate::to_global(momentum, normal);
}

Flux Flux::in_global(const Vector3d &normal) const {
    Flux f(*this);
    f.to_global(normal);
    return f;
}

std::ostream &operator<<(std::ostream &os, const Flux &flux) {
    os << "mass: " << flux.mass <<
       " momentum: {" << flux.momentum.x() << ", " << flux.momentum.y() << ", " << flux.momentum.z() <<
       "} energy: " << flux.energy;
    return os;
}

} // namespace smf


namespace mmf {

PState::PState()
    : density(0.0),
      velocity({0.0, 0.0, 0.0}),
      pressure(0.0),
      temperature(0.0),
      energy(0.0),
      mass_frac(),
      vol_frac() { }

PState::PState(double density, const Vector3d &velocity, double pressure,
        double energy, const Fractions &mass_frac,
        double temperature, const Fractions& vol_frac)
    : density(density),
      velocity(velocity),
      pressure(pressure),
      energy(energy),
      temperature(temperature),
      mass_frac(mass_frac),
      vol_frac(vol_frac) {
}

PState::PState(const QState &q, const phys::Materials &mixture, double P0, double T0) {
    mass_frac = Fractions(q.mass_frac);

    density = q.mass;
    velocity = q.momentum / density;
    energy = q.energy / density - 0.5 * velocity.squaredNorm();
    pressure = mixture.pressure_re(density, energy, mass_frac, {.P0=P0, .T0=T0});
    if (std::isnan(pressure)) {
        pressure = mixture.pressure_re(density, energy, mass_frac, {.P0=-P0, .T0=T0});
    }
    temperature = mixture.temperature_rp(density, pressure, mass_frac, {.T0=T0});
}

void PState::to_local(const Vector3d &normal) {
    Rotate::to_local(velocity, normal);
}

PState PState::in_local(const Vector3d &normal) const {
    PState z(*this);
    z.to_local(normal);
    return z;
}

void PState::to_global(const Vector3d &normal) {
    Rotate::to_global(velocity, normal);
}

PState PState::in_global(const Vector3d &normal) const {
    PState z(*this);
    z.to_global(normal);
    return z;
}

std::ostream &operator<<(std::ostream &os, const PState &state) {
    os << boost::format(
            "density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, temperature: %6%, energy: %7%, mass_frac: %8%") %
          state.density % state.velocity.x() % state.velocity.y() % state.velocity.z() % state.pressure %
          state.temperature % state.energy % state.mass_frac;
    return os;
}

ScalarSet PState::densities() const {
    ScalarSet densities;
    for (size_t i = 0; i < Fractions::size(); i++) {
        densities[i] = density * mass_frac[i];
    }

    return densities;
}

ScalarSet PState::energies() const {
    ScalarSet energies;
    for (size_t i = 0; i < Fractions::size(); i++) {
        energies[i] = energy * mass_frac[i];
    }

    return energies;
}

smf::PState PState::to_smf() const {
    return {density, velocity, pressure, energy};
}

bool PState::is_bad() const {
    return std::isinf(density) || std::isnan(density) ||
           std::isinf(velocity.x()) || std::isnan(velocity.x()) ||
           std::isinf(velocity.y()) || std::isnan(velocity.y()) ||
           std::isinf(velocity.z()) || std::isnan(velocity.z()) ||
           std::isinf(pressure) || std::isnan(pressure) ||
           std::isinf(energy) || std::isnan(energy) ||
           std::isinf(temperature) || std::isnan(temperature) ||
           mass_frac.empty() ||
           density < 0;
}


QState::QState()
    : mass(0.0),
      momentum({0.0, 0.0, 0.0}),
      energy(0.0),
      mass_frac() {

}

QState::QState(double mass, const Vector3d &momentum, double energy, const ScalarSet &mass_frac)
    : mass(mass),
      momentum(momentum),
      energy(energy),
      mass_frac(mass_frac) {

}

QState::QState(const PState &z) {
    mass_frac = ScalarSet(z.mass_frac);
    mass_frac *= z.density;

    mass = z.density;
    momentum = z.density * z.velocity;
    energy = z.density * (z.energy + 0.5 * z.velocity.squaredNorm());
}

void QState::to_local(const Vector3d &normal) {
    Rotate::to_local(momentum, normal);
}

QState QState::in_local(const Vector3d &normal) const {
    QState q(*this);
    q.to_local(normal);
    return q;
}

void QState::to_global(const Vector3d &normal) {
    Rotate::to_global(momentum, normal);
}

QState QState::in_global(const Vector3d &normal) const {
    QState q(*this);
    q.to_global(normal);
    return q;
}

std::ostream &operator<<(std::ostream &os, const QState &state) {
    os << boost::format("mass: %1%, momentum: {%2%, %3%, %4%}, energy: %5%, mass_frac: %6%") %
          state.mass %
          state.momentum.x() % state.momentum.y() % state.momentum.z() %
          state.energy %
          state.mass_frac;
    return os;
}

Flux::Flux()
    : mass(0.0),
      momentum({0.0, 0.0, 0.0}),
      energy(0.0),
      mass_frac() {

}

Flux::Flux(const PState &z) {
    mass_frac = ScalarSet(z.mass_frac);
    mass_frac *= z.density * z.velocity.x();

    mass = z.density * z.velocity.x();
    momentum.x() = z.density * z.velocity.x() * z.velocity.x() + z.pressure;
    momentum.y() = z.density * z.velocity.x() * z.velocity.y();
    momentum.z() = z.density * z.velocity.x() * z.velocity.z();
    energy = (z.density * (z.energy + 0.5 * z.velocity.squaredNorm()) + z.pressure) * z.velocity.x();
}

void Flux::to_local(const Vector3d &normal) {
    Rotate::to_local(momentum, normal);
}

Flux Flux::in_local(const Vector3d &normal) const {
    Flux f(*this);
    f.to_local(normal);
    return f;
}

void Flux::to_global(const Vector3d &normal) {
    Rotate::to_global(momentum, normal);
}

Flux Flux::in_global(const Vector3d &normal) const {
    Flux f(*this);
    f.to_global(normal);
    return f;
}

std::ostream &operator<<(std::ostream &os, const Flux &flux) {
    os << boost::format("mass: %1%, momentum: {%2%, %3%, %4%}, energy: %5%, mass_frac: %6%") %
          flux.mass %
          flux.momentum.x() % flux.momentum.y() % flux.momentum.z() %
          flux.energy %
          flux.mass_frac;
    return os;
}

Flux::Flux(double mass, const Vector3d &momentum,
        double energy, const ScalarSet &mass_frac)
    : mass(mass),
      momentum(momentum),
      energy(energy),
      mass_frac(mass_frac) {

}

bool Flux::is_bad() const {
    return std::isinf(mass) || std::isnan(mass) ||
           std::isinf(momentum.x()) || std::isnan(momentum.x()) ||
           std::isinf(momentum.y()) || std::isnan(momentum.y()) ||
           std::isinf(momentum.z()) || std::isnan(momentum.z()) ||
           std::isinf(energy) || std::isnan(energy);
}

} // namespace mmf

} // namespace zephyr