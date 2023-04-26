#include <zephyr/math/cfd/models.h>

namespace zephyr {
namespace math {

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
    os << "density: " << state.density << " velocity: {" << state.velocity.x() << ", " << state.velocity.y() << ", " << state.velocity.z() <<
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

Flux::Flux()
        : mass(0.0), momentum(0.0, 0.0, 0.0), energy(0.0) {

}

Flux::Flux(const PState &state) {
    mass = state.density * state.velocity.x();
    momentum.x() = state.density * state.velocity.x() * state.velocity.x() + state.pressure;
    momentum.y() = state.density * state.velocity.x() * state.velocity.y();
    momentum.z() = state.density * state.velocity.x() * state.velocity.z();
    energy = (state.density * (state.energy + 0.5 * state.velocity.squaredNorm()) + state.pressure) * state.velocity.x();
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
    os << "mass: " << flux.mass << " momentum: {" << flux.momentum.x() << ", " << flux.momentum.y() << ", " << flux.momentum.z() <<
       "} energy: " << flux.energy;
    return os;
}

Flux::Flux(double mass, const Vector3d &momentum, double energy) : mass(mass), momentum(momentum), energy(energy) {}

} // namespace smf

} // namespace math
} // namespace zephyr