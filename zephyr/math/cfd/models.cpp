#include <zephyr/math/cfd/models.h>
#include <boost/format.hpp>

namespace zephyr::math {

namespace smf {

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
        : mass(mass), momentum(momentum), energy(energy) {
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

Flux::Flux() : mass(0.0), momentum(0.0, 0.0, 0.0), energy(0.0) {}

Flux::Flux(const PState &state) {
    mass = state.density * state.velocity.x();
    momentum.x() = state.density * state.velocity.x() * state.velocity.x() + state.pressure;
    momentum.y() = state.density * state.velocity.x() * state.velocity.y();
    momentum.z() = state.density * state.velocity.x() * state.velocity.z();
    energy =
            (state.density * (state.energy + 0.5 * state.velocity.squaredNorm()) + state.pressure) * state.velocity.x();
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

Flux::Flux(double mass, const Vector3d &momentum, double energy) : mass(mass), momentum(momentum), energy(energy) {}

} // namespace smf


namespace mmf {

PState::PState(const double &pressure, const double &temperature,
               const Vector3d &velocity, const std::vector<Component> &components) :
        pressure(pressure), temperature(temperature), velocity(velocity) {
    if (components.size() > Fractions::max_size) {
        throw std::runtime_error("When construct PState got components.size() > Fractions::max_size (" +
                                 std::to_string(components.size()) + " > " + std::to_string(Fractions::max_size));
    }

    std::vector<double> fracs(components.size());
    density = 0;
    for (size_t i = 0; i < components.size(); ++i) {
        fracs[i] = components[i].frac;
    }
    mass_frac = Fractions(fracs);

    energy = 0;
    for (size_t i = 0; i < components.size(); ++i) {
        density += mass_frac[i] * components[i].density;
        energy += mass_frac[i] * components[i].energy;
    }
}

PState::PState(const double &density, const Vector3d &velocity,
               const double &pressure, const double &energy, const double &temperature, const Fractions &mass_frac)
        : density(density), velocity(velocity),
          pressure(pressure), energy(energy), temperature(temperature), mass_frac(mass_frac) {
}

PState::PState(const QState &q, const phys::Materials &mixture, double P0, double T0) {
    density = q.mass;
    velocity = q.momentum / density;
    energy = q.energy / density - 0.5 * velocity.squaredNorm();

    mass_frac = Fractions(q.mass_frac);

    pressure = mixture.pressure_re(density, energy, mass_frac, {.P0=P0});
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

std::vector<double> PState::get_densities() const {
    std::vector<double> densities(mass_frac.get_size());
    for (size_t i = 0; i < mass_frac.get_size(); i++) {
        densities[i] = density * mass_frac[i];
    }

    return densities;
}

std::vector<double> PState::get_energies() const {
    std::vector<double> energies(mass_frac.get_size());
    for (size_t i = 0; i < mass_frac.get_size(); i++) {
        energies[i] = density * mass_frac[i];
    }

    return energies;
}

smf::PState PState::to_smf() const {
    return {density, velocity, pressure, energy};
}


QState::QState(const double &mass, const Vector3d &momentum, const double &energy, const FractionsFlux &mass_frac)
        : mass(mass), momentum(momentum), energy(energy), mass_frac(mass_frac) {}

QState::QState(const PState &state) {
    mass_frac = FractionsFlux(state.mass_frac);
    mass_frac *= state.density;

    mass = state.density;
    momentum = state.density * state.velocity;
    energy = state.density * (state.energy + 0.5 * state.velocity.squaredNorm());
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

Flux::Flux() : mass(0.0), momentum(0.0, 0.0, 0.0), energy(0.0), mass_frac() {}

Flux::Flux(const PState &state) {
    mass_frac = FractionsFlux(state.mass_frac);
    mass_frac *= state.density * state.velocity.x();

    mass = state.density * state.velocity.x();
    momentum.x() = state.density * state.velocity.x() * state.velocity.x() + state.pressure;
    momentum.y() = state.density * state.velocity.x() * state.velocity.y();
    momentum.z() = state.density * state.velocity.x() * state.velocity.z();
    energy =
            (state.density * (state.energy + 0.5 * state.velocity.squaredNorm()) + state.pressure) * state.velocity.x();
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

Flux::Flux(double mass, const Vector3d &momentum, double energy,
           const FractionsFlux &mass_frac) : mass(mass), momentum(momentum), energy(energy), mass_frac(mass_frac) {}


} // namespace mmf

} // namespace zephyr