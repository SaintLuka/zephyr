#include <iostream>
#include <algorithm>

#include <zephyr/geom/vector.h>
#include <zephyr/math/cfd/flux/hll.h>

namespace zephyr::math {

smf::Flux HLL::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux HLL::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    // Нормальные скорости слева и справа
    const double& u_L = zL.velocity.x();
    const double& u_R = zR.velocity.x();

    // Скорость звука слева и справа
    double c_L = eos.sound_speed_rP(zL.density, zL.pressure);
    double c_R = eos.sound_speed_rP(zR.density, zR.pressure);

    // Оценки скоростей расходящихся волн
    double S_L = std::min({u_L - c_L, u_R - c_R, 0.0});
    double S_R = std::max({u_L + c_L, u_R + c_R, 0.0});

    QState Q_L(zL); // Консервативный вектор слева
    QState Q_R(zR); // Консервативный вектор справа

    Flux F_L(zL);   // Дифференциальный поток слева
    Flux F_R(zR);   // Дифференциальный поток справа

    Flux F = (S_R * F_L.arr() - S_L * F_R.arr() + S_L * S_R * (Q_R.arr() - Q_L.arr())) / (S_R - S_L);

    if (F.arr().hasNaN()) {
        std::cerr << "HLL::calc_flux error\n";
        std::cerr << "  z_L: " << zL << "\n";
        std::cerr << "  z_R: " << zR << "\n";
        std::cerr << "  c_L: " << c_L << "; c_R: " << c_R << "\n";
        std::cerr << "  S_L: " << S_L << "; S_R: " << S_R << "\n";
        std::cerr << "  F_HLL: " << F << "\n";
        throw std::runtime_error("HLL::calc_flux error: bad value");
    }

    return F;
}

HLL::WaveConfig HLL::wave_config(const phys::Eos& eos,
        const smf::PState& zL, const smf::PState& zR) {
    // Нормальные скорости слева и справа
    const double &u_L = zL.velocity.x();
    const double &u_R = zR.velocity.x();

    // Скорость звука слева и справа
    double c_L = eos.sound_speed_rP(zL.density, zL.pressure);
    double c_R = eos.sound_speed_rP(zR.density, zR.pressure);

    // Оценки скоростей расходящихся волн
    double S_L = std::min({u_L - c_L, u_R - c_R, 0.0});
    double S_R = std::max({u_L + c_L, u_R + c_R, 0.0});

    smf::QState Q_L(zL); // Консервативный вектор слева
    smf::QState Q_R(zR); // Консервативный вектор справа

    smf::Flux F_L(zL);   // Дифференциальный поток слева
    smf::Flux F_R(zR);   // Дифференциальный поток справа

    smf::Flux F = (S_R * F_L.arr() - S_L * F_R.arr() + S_L * S_R * (Q_R.arr() - Q_L.arr())) / (S_R - S_L);
    smf::QState Q = (F_L.arr() - F_R.arr() + S_R * Q_R.arr() - S_L * Q_L.arr()) / (S_R - S_L);

    return WaveConfig({.S_L = S_L, .S_R = S_R, .Qs = Q, .Fs = F});
}

HLL::WaveConfig HLL::wave_config(const phys::Eos& eos,
                            const smf::QState& Q_L, const smf::Flux& F_L,
                            const smf::QState& Q_R, const smf::Flux& F_R) {
    Vector3d vL = Q_L.momentum / Q_L.density;
    Vector3d vR = Q_R.momentum / Q_R.density;

    double eL = Q_L.energy / Q_L.density - 0.5 * vL.squaredNorm();
    double eR = Q_R.energy / Q_R.density - 0.5 * vR.squaredNorm();

    // Скорость звука слева и справа
    double c_L = eos.sound_speed_re(Q_L.density, eL);
    double c_R = eos.sound_speed_re(Q_R.density, eR);

    // Оценки скоростей расходящихся волн
    double S_L = std::min({vL.x() - c_L, vR.x() - c_R, 0.0});
    double S_R = std::max({vL.x() + c_L, vR.x() + c_R, 0.0});

    // Поток через левую/правую волну
    smf::Flux a_L = F_L.arr() - S_L * Q_L.arr();
    smf::Flux a_R = F_R.arr() - S_R * Q_R.arr();

    smf::Flux F = (S_R * a_L.arr() - S_L * a_R.arr()) / (S_R - S_L);
    smf::QState Q = (a_L.arr() - a_R.arr()) / (S_R - S_L);

#if 0
    double err1 = ((F_R.arr() - F.arr()) - S_R * (Q_R.arr() - Q.arr())).cwiseAbs().maxCoeff();
    double err2 = ((F.arr() - F_L.arr()) - S_L * (Q.arr() - Q_L.arr())).cwiseAbs().maxCoeff();

    if (std::max({err1, err2}) > 1.0e-6) {
        throw std::runtime_error("Rankine-Hugoniot conditions failed");
    }
#endif

    return WaveConfig({.S_L = S_L, .S_R = S_R, .Qs = Q, .Fs = F});
}

mmf::Flux HLL::flux(const mmf::PState &zL, const mmf::PState &zR, const phys::MixturePT &mixture) const {
    return calc_flux(zL, zR, mixture);
}

mmf::Flux HLL::calc_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::MixturePT &mixture) {
    using namespace mmf;

    // Нормальные скорости слева и справа
    const double& u_L = zL.velocity.x();
    const double& u_R = zR.velocity.x();

    // Скорость звука слева и справа
    double c_L = mixture.sound_speed_rP(zL.density, zL.pressure, zL.mass_frac,
                                        {.T0 = zL.T(), .rhos = &zL.densities});
    double c_R = mixture.sound_speed_rP(zR.density, zR.pressure, zR.mass_frac,
                                        {.T0 = zR.T(), .rhos = &zR.densities});

    // Оценки скоростей расходящихся волн
    double S_L = std::min({u_L - c_L, u_R - c_R, 0.0});
    double S_R = std::max({u_L + c_L, u_R + c_R, 0.0});

    QState Q_L(zL); // Консервативный вектор слева
    QState Q_R(zR); // Консервативный вектор справа

    Flux F_L(zL);   // Дифференциальный поток слева
    Flux F_R(zR);   // Дифференциальный поток справа

    Flux F = (S_R * F_L.arr() - S_L * F_R.arr() + S_L * S_R * (Q_R.arr() - Q_L.arr())) / (S_R - S_L);

    if (F.arr().hasNaN()) {
        std::cerr << "HLL::calc_flux error\n";
        std::cerr << "  z_L: " << zL << "\n";
        std::cerr << "  z_R: " << zR << "\n";
        std::cerr << "  c_L: " << c_L << "; c_R: " << c_R << "\n";
        std::cerr << "  S_L: " << S_L << "; S_R: " << S_R << "\n";
        std::cerr << "  F_HLL: " << F << "\n";
        throw std::runtime_error("HLL::calc_flux error: bad value");
    }

    return F;
}

} // namespace zephyr::math