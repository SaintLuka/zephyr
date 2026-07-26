#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <format>
#include <boost/format.hpp>
#include <zephyr/math/funcs.h>
#include <zephyr/math/cfd/models.h>

namespace zephyr::math {

namespace swe {

PState::PState() : depth(0.0), velocity({0.0, 0.0}) {}

PState::PState(double depth, const Vector2d &velocity)
    : depth(depth), velocity(velocity) {
    z_assert(depth >= 0.0);
}

PState::PState(const QState &q) {
    if (q.depth <= 0.0) {
        depth = 0.0;
        velocity = {0.0, 0.0};
    }
    else {
        depth = q.depth;
        velocity = q.momentum / q.depth;
    }
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

bool PState::is_bad() const {
    return !std::isfinite(depth) || depth < 0.0 || !std::isfinite(velocity.x()) || !std::isfinite(velocity.y());
}

std::ostream &operator<<(std::ostream &os, const PState &state) {
    os << std::format("h: {:.5f},  v: [{:.5e}, {:.5e}]",
        state.depth, state.velocity.x(), state.velocity.y());
    return os;
}

QState::QState() : depth(0.0), momentum({0.0, 0.0}) {}

QState::QState(double depth, const Vector2d &momentum)
    : depth(depth), momentum(momentum) {}

QState::QState(const PState &z)
    : depth(z.depth), momentum(z.velocity * depth) {
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
    os << std::format("h: {:.5f},  hv: [{:.5e}, {:.5e}]",
        state.depth, state.momentum.x(), state.momentum.y());
    return os;
}

Flux::Flux() : mass(0.0), momentum({0.0, 0.0}) { }

Flux::Flux(double mass, const Vector2d &momentum)
    : mass(mass), momentum(momentum) { }

Flux::Flux(const PState& z) {
    mass = z.depth * z.velocity.x();
    momentum.x() = mass * z.velocity.x() + 0.5 * g * z.depth * z.depth;
    momentum.y() = mass * z.velocity.y();
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
    os << std::format("h: {:.5f},  hv: [{:.5e}, {:.5e}]",
        flux.mass, flux.momentum.x(), flux.momentum.y());
    return os;
}

} // namespace swe

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

bool PState::is_bad(const phys::Eos &eos) const {
    if (!std::isfinite(density) || !std::isfinite(pressure) || !std::isfinite(energy)) {
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

    std::tie(densities, energy, temperature)
        = mixture.get_reT(density, pressure, mass_frac);
}

PState::PState(const QState &q, const phys::MixturePT &mixture,
               double P0, double T0, std::span<const double> rhos0) {

    density   = q.density;
    velocity  = q.momentum / density;
    energy    = q.energy / density - 0.5 * velocity.squaredNorm();
    for (int i = 0; i < mixture.size(); ++i) {
        mass_frac[i] = q.mass_frac[i] / density;
    }
    mass_frac.cutoff(1.0e-12);
    mass_frac.normalize();

    std::tie(densities, pressure, temperature) =
        mixture.get_rPT(density, energy, mass_frac, {.P0=P0, .T0=T0, .rhos=rhos0});
}

PState PState::Mix1(
        const MixturePT &mixture,
        const std::vector<double>& vol_fracs,
        const std::vector<PState> &zs,
        std::vector<int> indices) {

    if (vol_fracs.size() < 2) {
        throw std::runtime_error("PState::Mix1: require at least two volume fractions");
    }
    if (vol_fracs.size() != zs.size()) {
        throw std::runtime_error("PState::Mix1: number of volume fractions != number of states");
    }

    int n = static_cast<int>(vol_fracs.size());

    if (indices.empty()) {
        if (n > mixture.size()) {
            throw std::runtime_error("PState::Mix1: number of volume fractions more than materials");
        }
        indices.resize(n);
        for (int i = 0; i < n; ++i) {
            indices[i] = i;
        }
    }
    else {
        if (indices.size() != n) {
            throw std::runtime_error("PState::Mix1: number of volume fractions != number of indices");
        }
        for (int id: indices) {
            if (id >= mixture.size()) {
                throw std::runtime_error("PState::Mix1: material index >= mixture.size()");
            }
        }
    }

    // Копируем и нормируем массив объемных долей
    Fractions alpha = Fractions::Zero();
    for (int i = 0; i < n; ++i) {
        alpha[indices[i]] += vol_fracs[i];
    }
    alpha.normalize();

    double temperature = 0.0;
    for (int i = 0; i < n; ++i) {
        // Вариант 1. Усреднение температуры по объемным долям
        temperature += alpha[indices[i]] * zs[i].temperature;

        // Вариант 2. Минимум по всем состояниям
        // temperature = std::min(temperature, zs[i].temperature);
    }

    // Проверяем совпадение всех давлений
    double pressure = zs[0].pressure;
    for (int i = 1; i < n; ++i) {
        if (std::abs(zs[i].pressure - pressure) > 1.0e-10 * std::abs(pressure)) {
            throw std::runtime_error("PState::Mix1: assuming pressures are the same for each state");
        }
    }

    // Истинные плотности компонент (отличаются от исходных)
    ScalarSet densities = ScalarSet::NaN();
    for (int i = 0; i < mixture.size(); ++i) {
        if (alpha.has(i)) {
            densities[i] = 1.0 / mixture[i].volume_PT(pressure, temperature);
        }
    }

    // Средняя плотность
    double density = 0.0;
    for (int i = 0; i < mixture.size(); ++i) {
        if (alpha.has(i)) {
            density += alpha[i] * densities[i];
        }
    }

    // Массовые доли
    Fractions mass_frac = Fractions::Zero();
    for (int i = 0; i < mixture.size(); ++i) {
        if (alpha.has(i)) {
            mass_frac[i] = alpha[i] * densities[i] / density;
        }
    }
    mass_frac.normalize();

    // Взвешенная скорость и внутренняя энергия
    Vector3d velocity = Vector3d::Zero();
    for (int i = 0; i < n; ++i) {
        if (alpha.has(indices[i])) {
            velocity += mass_frac[indices[i]] * zs[i].velocity;
        }
    }

    return PState(density, velocity, pressure, mass_frac, mixture);
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

double PState::alpha(int i) const {
    return std::isnan(densities[i]) ? 0.0 :
           std::max(0.0, std::min(mass_frac[i] * density / densities[i], 1.0));
}

Fractions PState::volume_fractions() const {
    Fractions alpha = Fractions::Zero();
    for (int i = 0; i < mass_frac.size(); ++i) {
        if (mass_frac.has(i)) {
            // Если есть массовая концентрация beta_i, значит
            // должна быть определена плотность rho_i !
            alpha[i] = mass_frac[i] * density / densities[i];
        }
    }
    return alpha;
}

std::ostream &operator<<(std::ostream &os, const PState &state) {
    os << boost::format(
            "ρ: %.5f,  v: {%+.5e, %+.5e, %+.5e},  P: %+.5e,  e: %+.5e,  T: %+.5e,  ") %
          state.density % state.velocity.x() % state.velocity.y() % state.velocity.z() %
          state.pressure % state.energy % state.temperature;
    os << boost::format("β: %1%,  ϱ: %2%") % state.mass_frac % state.densities;
    return os;
}

double PState::true_energy(const MixturePT& mixture, int i) const {
    return mixture[i].energy_rT(densities[i], temperature, {.deriv = false});
}

smf::PState PState::to_smf() const {
    return {density, velocity, pressure, energy};
}

smf::PState PState::extract(const MixturePT& mixture, int idx) const {
    return smf::PState(densities[idx], velocity, pressure, true_energy(mixture, idx));
}

std::pair<mmf::PState, mmf::PState> PState::split(const MixturePT& mixture, int iA) const {
    // Внутренняя энергия для материала A
    double energy_A = mixture[iA].energy_rT(densities[iA], temperature, {.P0=pressure});

    // Чистое состояние для материала A
    mmf::PState zA(
            densities[iA],
            velocity,
            pressure,
            energy_A,
            temperature,
            Fractions::Pure(iA),
            ScalarSet::PureNaN(iA, densities[iA]));

    // Смешанное состояние (всё кроме A)
    mmf::PState zB(
            NAN,
            velocity,
            pressure,
            NAN,
            temperature,
            Fractions::Zero(),
            ScalarSet::NaN());

#if 1 // MRV VERSION

    double beta_A = mass_frac[iA];
    double beta_B = 1.0 - mass_frac[iA];

    // Для работы функции обязательно выполнение условия
    // 1/rho = sum_i beta_i / rho_i
    // Если нет совместности, то split будет приводить к ошибкам.
    // zB.density = beta_B / (1.0 / density - beta_A / densities[iA]);

    // Эта версия работает, даже если не выполнено условие совместности.
    double denom = 0.0;
    for (int i = 0; i < mass_frac.size(); ++i) {
        if ( mass_frac.has(i) && i != iA) {
            denom += mass_frac[i] / densities[i];
        }
    }
    zB.density = beta_B / denom;

    zB.energy  = (energy - beta_A * energy_A ) / beta_B;

    for (int i = 0; i < mass_frac.size(); ++i) {
        if (mass_frac.has(i) && i != iA) {
            zB.densities[i] = densities[i];
            zB.mass_frac[i] = mass_frac[i] / beta_B;
        }
    }

    if (zB.is_bad()) {
        std::cout << "bad split #1\n";
    }

#else

    double beta_B = 1.0 - mass_frac[iA];
    Fractions alpha;
    for (int i = 0; i < mass_frac.size(); ++i) {
        if( mass_frac.has(i) && i != iA ) {
            zB.densities[i] = densities[i];
            zB.mass_frac[i] = mass_frac[i] / beta_B;
        }
        if( mass_frac.has(i) ) {
            alpha[i] = density*mass_frac[i]/densities[i];
            alpha[i] = this->alpha(i);
        }
    }
    zB.mass_frac.normalize();

    //alphas at splitted z2
    double alpha_B = 1.0 - alpha[iA];
    alpha[iA] = 0.0;
    for (int i = 0; i < mass_frac.size(); ++i) {
        if( mass_frac.has(i) && i != iA ) {
            alpha[i] = alpha[i]/alpha_B;
        }
    }
    alpha.normalize();
    //mixture density at splited zB
    zB.density = 0.0;
    for (int i = 0; i < mass_frac.size(); ++i) {
        if( mass_frac.has(i) && i != iA ) {
            zB.density = zB.density + densities[i]*alpha[i];
        }
    }
    //betas at splitted zB
    zB.energy = 0.0;
    for (int i = 0; i < mass_frac.size(); ++i) {
        if( mass_frac.has(i) && i != iA ) {
            zB.mass_frac[i]  = alpha[i]*densities[i]/zB.density;
            zB.energy += zB.mass_frac[i]*mixture[i].energy_rT(densities[i], temperature, {.P0=pressure});
        }
    }
    if(zB.is_bad()) {
        std::cout << "bad split #2\n";
    }

#endif
    return {zA, zB};
}

void PState::interpolation_update(const MixturePT& mixture) {
    // Нормализуем после интерполяции
    mass_frac.normalize();

    // Восстанавливаем совместность после интерполяции
    auto[rhos, e, T] = mixture.get_reT(density, pressure, mass_frac,
                                       {.T0=temperature, .rhos=densities.span()});
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
    density   = z.rho();
    momentum  = z.rho() * z.vel();
    energy    = z.rho() * z.E();
    for (int i = 0; i < mass_frac.size(); ++i) {
        mass_frac[i] = z.rho() * z.beta(i);
    }
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
    density      = z.rho() * z.vx();
    momentum.x() = z.rho() * z.vx() * z.vx() + z.P();
    momentum.y() = z.rho() * z.vx() * z.vy();
    momentum.z() = z.rho() * z.vx() * z.vz();
    energy = (z.rho() * z.E() + z.P()) * z.vx();
    for (int i = 0; i < mass_frac.size(); ++i) {
        mass_frac[i] = z.rho() * z.vx() * z.beta(i);
    }
}

Flux::Flux(const smf::Flux& flux, int mat)
    : density(flux.density),
      momentum(flux.momentum),
      energy(flux.energy),
      mass_frac(ScalarSet::Pure(mat, flux.density)) {
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
    for (int i = 0; i < mass_frac.size(); ++i) {
        mass_frac[i] = -mass_frac[i];
    }
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