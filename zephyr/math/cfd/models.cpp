#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <boost/format.hpp>

#include <zephyr/math/cfd/models.h>

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
    density = q.density;
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

void PState::inverse() {
    velocity.x() = -velocity.x();
}

bool PState::is_bad(const phys::Eos &eos) {
    if (std::isnan(density) || std::isnan(pressure) || std::isnan(energy)) {
        return true;
    }
    return density <= 0.0 || pressure < eos.min_pressure();
}

std::ostream &operator<<(std::ostream &os, const PState &state) {
    os << boost::format("ρ: %.5f,  v: {%+.5e, %+.5e, %+.5e},  P: %+.5e,  e: %+.5e") %
          state.density % state.velocity.x() % state.velocity.y() % state.velocity.z() %
          state.pressure % state.energy;
    return os;
}

QState::QState()
    : density(0.0),
      momentum({0.0, 0.0, 0.0}),
      energy(0.0) {

}

QState::QState(const double &mass, const Vector3d &momentum, const double &energy)
    : density(mass),
      momentum(momentum),
      energy(energy) {
}

QState::QState(const PState &z) {
    density = z.density;
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
    os << boost::format("ρ: %.5f,  ρv: {%+.5e, %+.5e, %+.5e},  ρE: %+.5e") %
          state.density % state.momentum.x() % state.momentum.y() % state.momentum.z() % state.energy;
    return os;
}

Flux::Flux()
    : density(0.0),
      momentum({0.0, 0.0, 0.0}),
      energy(0.0) {

}

Flux::Flux(double mass, const Vector3d &momentum, double energy)
    : density(mass),
      momentum(momentum),
      energy(energy) {

}

Flux::Flux(const PState &z) {
    density = z.density * z.velocity.x();
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
    os << boost::format("ρ: %+.5f,  ρv: {%+.5e, %+.5e, %+.5e},  ρE: %+.5e") %
          flux.density % flux.momentum.x() % flux.momentum.y() % flux.momentum.z() % flux.energy;
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
      mass_frac{0.0},
      densities{0.0} { }

PState::PState(double density, const Vector3d &velocity, double pressure,
        double energy, double temperature,
        const Fractions &mass_frac, const ScalarSet& rhos)
    : density(density),
      velocity(velocity),
      pressure(pressure),
      energy(energy),
      temperature(temperature),
      mass_frac(mass_frac),
      densities(rhos) {
}

PState::PState(double density, const Vector3d &velocity, double pressure,
       const Fractions &mass_frac,  const MixturePT &mixture)
    : density(density),
      velocity(velocity),
      pressure(pressure),
      mass_frac(mass_frac) {

    auto[rhos, e, T] = mixture.get_reT(density, pressure, mass_frac);
    densities   = rhos;
    energy      = e;
    temperature = T;
}

PState::PState(const QState &q, const phys::MixturePT &mixture,
               double P0, double T0, const ScalarSet& rhos0) {

    density   = q.density;
    velocity  = q.momentum / density;
    energy    = q.energy / density - 0.5 * velocity.squaredNorm();

    mass_frac = q.mass_frac.arr() / density;
    mass_frac.cutoff(1.0e-12);
    mass_frac.normalize();

    auto[rhos, P, T] = mixture.get_rPT(density, energy, mass_frac,
                                       {.P0=P0, .T0=T0, .rhos=&rhos0});

    densities   = rhos;
    pressure    = P;
    temperature = T;
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

void PState::inverse() {
    velocity.x() = -velocity.x();
}

std::ostream &operator<<(std::ostream &os, const PState &state) {
    os << boost::format(
            "ρ: %.5f,  v: {%+.5e, %+.5e, %+.5e},  P: %+.5e,  e: %+.5e,  T: %+.5e,  ") %
          state.density % state.velocity.x() % state.velocity.y() % state.velocity.z() %
          state.pressure % state.energy % state.temperature;
    os << boost::format("β: %1%,  ϱ: %2%") % state.mass_frac % state.densities;
    return os;
}

double PState::true_energy(const MixturePT& mixture, int idx) const {
    return mixture[idx].energy_rT(densities[idx], temperature, {.deriv = false});
}

smf::PState PState::to_smf() const {
    return {density, velocity, pressure, energy};
}

smf::PState PState::extract(const MixturePT& mixture, int idx) const {
    return smf::PState(densities[idx], velocity, pressure, true_energy(mixture, idx));
}

std::pair<mmf::PState, mmf::PState> PState::split(const MixturePT& mixture, int iA) const {
    // Внутренняя энегия для материала A
    double energy_A = mixture[iA].energy_rT(densities[iA], temperature, {.P0=pressure});

    // Чистое состояние для материала A
    mmf::PState zA(
            densities[iA],
            velocity,
            pressure,
            energy_A,
            temperature,
            Fractions::Pure(iA),
            ScalarSet::Pure(iA, densities[iA]));

    // Смешаное состояние (всё кроме A)
    mmf::PState zB(
            NAN,
            velocity,
            pressure,
            NAN,
            temperature,
            Fractions::Zero(),
            ScalarSet::NaN());

    // Массовые концентрации, с которыми смешиваются zA и zB
    double beta_A = mass_frac[iA];
    double beta_B = 1.0 - mass_frac[iA];

    zB.density = beta_B / (1.0 / density - beta_A / densities[iA]);
    zB.energy  = (energy - beta_A * energy_A ) / beta_B;

    for (int i = 0; i < Fractions::size(); ++i) {
        if (mass_frac.has(i) && i != iA) {
            zB.densities[i] = densities[i];
            zB.mass_frac[i] = mass_frac[i] / beta_B;
        }
    }

    // Случай небольших объемных долей. Формулы выше математически верны,
    // но работают неточно, что приводит к ошибкам.
    if (beta_B < 1.0e-12) {
        zB.mass_frac.normalize();

        double mix_vol = 0.0;
        double mix_e = 0.0;
        for (int i = 0; i < Fractions::size(); ++i) {
            if (zB.mass_frac.has(i)) {
                mix_vol += zB.mass_frac[i] / zB.densities[i];
                mix_e += zB.mass_frac[i] * mixture[i].energy_rT(zB.densities[i], zB.temperature);
            }
        }
        zB.density = 1.0 / mix_vol;
        zB.energy  = mix_e;
    }

    return {zA, zB};
}

void PState::interpolation_update(const MixturePT& mixture) {
    // Нормализуем после интерполяции
    mass_frac.normalize();

    // Восстанавливаем совместность после интерполяции
    auto[rhos, e, T] = mixture.get_reT(density, pressure, mass_frac,
                                       {.T0=temperature, .rhos=&densities});
    energy      = e;
    temperature = T;
    densities   = rhos;
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
    : density(0.0),
      momentum({0.0, 0.0, 0.0}),
      energy(0.0),
      mass_frac() {

}

QState::QState(double mass, const Vector3d &momentum, double energy, const ScalarSet &mass_frac)
    : density(mass),
      momentum(momentum),
      energy(energy),
      mass_frac(mass_frac) {

}

QState::QState(const PState &z) {
    density      = z.density;
    momentum  = z.density * z.velocity;
    energy    = z.density * (z.energy + 0.5 * z.velocity.squaredNorm());
    mass_frac = z.density * z.mass_frac.arr();
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
    os << boost::format("ρ: %.5f,  ρv: {%+.5e, %+.5e, %+.5e},  ρE: %+.5e,  ") %
          state.density %
          state.momentum.x() % state.momentum.y() % state.momentum.z() %
          state.energy;
    os << boost::format("ρβ: %1%") % state.mass_frac;
    return os;
}

Flux::Flux()
    : density(0.0),
      momentum({0.0, 0.0, 0.0}),
      energy(0.0),
      mass_frac() {

}

Flux::Flux(const PState &z) {
    density      = z.density * z.velocity.x();
    mass_frac    = z.density * z.velocity.x() * z.mass_frac.arr();
    momentum.x() = z.density * z.velocity.x() * z.velocity.x() + z.pressure;
    momentum.y() = z.density * z.velocity.x() * z.velocity.y();
    momentum.z() = z.density * z.velocity.x() * z.velocity.z();
    energy = (z.density * (z.energy + 0.5 * z.velocity.squaredNorm()) + z.pressure) * z.velocity.x();
}

Flux::Flux(const smf::Flux& flux, int mat)
    : density(flux.density),
      momentum(flux.momentum),
      energy(flux.energy),
      mass_frac(mat, flux.density) {
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

void Flux::inverse() {
    density = -density;
    //momentum.x() = -momentum.x();
    momentum.y() = -momentum.y();
    momentum.z() = -momentum.z();
    energy = -energy;
    mass_frac *= -1.0;
}

std::ostream &operator<<(std::ostream &os, const Flux &flux) {
    os << boost::format("ρ: %+.5f,  ρv: {%+.5e, %+.5e, %+.5e},  ρE: %+.5e,  ") %
          flux.density %
          flux.momentum.x() % flux.momentum.y() % flux.momentum.z() %
          flux.energy;
    os << boost::format("ρβ: %1%") % flux.mass_frac;
    return os;
}

Flux::Flux(double mass, const Vector3d &momentum,
        double energy, const ScalarSet &mass_frac)
    : density(mass),
      momentum(momentum),
      energy(energy),
      mass_frac(mass_frac) {

}

bool Flux::is_bad() const {
    return std::isinf(density) || std::isnan(density) ||
           std::isinf(momentum.x()) || std::isnan(momentum.x()) ||
           std::isinf(momentum.y()) || std::isnan(momentum.y()) ||
           std::isinf(momentum.z()) || std::isnan(momentum.z()) ||
           std::isinf(energy) || std::isnan(energy);
}

} // namespace mmf

} // namespace zephyr