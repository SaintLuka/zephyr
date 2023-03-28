#include <zephyr/math/cfd/models.h>

namespace zephyr { namespace math {

PState::PState(const double &density, const Vector3d &velocity,
       const double &pressure, const double &energy)
        : density(density), velocity(velocity),
          pressure(pressure), energy(energy) {
}

PState::PState(const QState& q, const phys::Eos& eos) {
    density  = q.mass;
    velocity = q.momentum / density;
    energy   = q.energy / density - 0.5 * velocity.squaredNorm();
    pressure = eos.pressure_re(density, energy);
}

void PState::to_local(const Vector3d& normal) {
    Rotate::to_local(velocity, normal);
}

PState PState::in_local(const Vector3d& normal) const {
    PState z(*this);
    z.to_local(normal);
    return z;
}

void PState::to_global(const Vector3d& normal) {
    Rotate::to_global(velocity, normal);
}

PState PState::in_global(const Vector3d& normal) const {
    PState z(*this);
    z.to_global(normal);
    return z;
}

QState::QState(const double &mass, const Vector3d &momentum, const double &energy)
        : mass(mass), momentum(momentum), energy(energy) {
}

QState::QState(const PState& z) {
    mass     = z.density;
    momentum = z.density * z.velocity;
    energy   = z.density * (z.energy + 0.5 * z.velocity.squaredNorm());
}

void QState::to_local(const Vector3d& normal) {
    Rotate::to_local(momentum, normal);
}

QState QState::in_local(const Vector3d& normal) const {
    QState q(*this);
    q.to_local(normal);
    return q;
}

void QState::to_global(const Vector3d& normal) {
    Rotate::to_global(momentum, normal);
}

QState QState::in_global(const Vector3d& normal) const {
    QState q(*this);
    q.to_global(normal);
    return q;
}

Flux::Flux()
    : mass(0.0), momentum(0.0, 0.0, 0.0), energy(0.0) {

}

Flux::Flux(const PState& z) {
    mass = z.density * z.velocity.x();
    momentum.x() = z.density * z.velocity.x() * z.velocity.x() + z.pressure;
    momentum.y() = z.density * z.velocity.x() * z.velocity.y();
    momentum.z() = z.density * z.velocity.x() * z.velocity.z();
    energy = (z.density * (z.energy + 0.5 * z.velocity.squaredNorm()) + z.pressure) * z.velocity.x();
}

void Flux::to_local(const Vector3d& normal) {
    Rotate::to_local(momentum, normal);
}

Flux Flux::in_local(const Vector3d& normal) const {
    Flux f(*this);
    f.to_local(normal);
    return f;
}

void Flux::to_global(const Vector3d& normal) {
    Rotate::to_global(momentum, normal);
}

Flux Flux::in_global(const Vector3d& normal) const {
    Flux f(*this);
    f.to_global(normal);
    return f;
}

}
}