#include <iostream>

#include <zephyr/geom/vector.h>
#include <zephyr/math/cfd/flux/hll.h>
#include <zephyr/phys/fractions.h>

using namespace zephyr::phys;

namespace zephyr::math {

inline double min(double x, double y, double z) {
    return std::min(x, std::min(y, z));
}

inline double max(double x, double y, double z) {
    return std::max(x, std::max(y, z));
}

smf::Flux HLL::flux(const smf::PState &zL, const smf::PState &zR, const Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux HLL::calc_flux(const smf::PState &zL, const smf::PState &zR, const Eos &eos) {
    using namespace smf;

    // Нормальные скорости слева и справа
    const double& u_L = zL.velocity.x();
    const double& u_R = zR.velocity.x();

    // Скорость звука слева и справа
    double c_L = eos.sound_speed_rP(zL.density, zL.pressure);
    double c_R = eos.sound_speed_rP(zR.density, zR.pressure);

    // Оценки скоростей расходящихся волн
    double S_L = min(u_L - c_L, u_R - c_R, 0.0);
    double S_R = max(u_L + c_L, u_R + c_R, 0.0);

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

smf::WaveConfig2 HLL::wave_config(const Eos& eos,
        const smf::PState& zL, const smf::PState& zR) {
    // Нормальные скорости слева и справа
    const double &u_L = zL.velocity.x();
    const double &u_R = zR.velocity.x();

    // Скорость звука слева и справа
    double c_L = eos.sound_speed_rP(zL.density, zL.pressure);
    double c_R = eos.sound_speed_rP(zR.density, zR.pressure);

    // Оценки скоростей расходящихся волн
    double S_L = min(u_L - c_L, u_R - c_R, 0.0);
    double S_R = max(u_L + c_L, u_R + c_R, 0.0);

    smf::QState Q_L(zL); // Консервативный вектор слева
    smf::QState Q_R(zR); // Консервативный вектор справа

    smf::Flux F_L(zL);   // Дифференциальный поток слева
    smf::Flux F_R(zR);   // Дифференциальный поток справа

    smf::Flux F = (S_R * F_L.arr() - S_L * F_R.arr() + S_L * S_R * (Q_R.arr() - Q_L.arr())) / (S_R - S_L);
    smf::QState Q = (F_L.arr() - F_R.arr() + S_R * Q_R.arr() - S_L * Q_L.arr()) / (S_R - S_L);

    return smf::WaveConfig2({.S_L = S_L, .S_R = S_R, .Qs = Q, .Fs = F});
}

smf::WaveConfig2 HLL::wave_config(const Eos& eos,
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
    double S_L = min(vL.x() - c_L, vR.x() - c_R, 0.0);
    double S_R = max(vL.x() + c_L, vR.x() + c_R, 0.0);

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

    return smf::WaveConfig2({.S_L = S_L, .S_R = S_R, .Qs = Q, .Fs = F});
}

mmf::Flux HLL::flux(const mmf::PState &zL, const mmf::PState &zR, const MixturePT &mix) const {
    return calc_flux(zL, zR, mix);
}

mmf::Flux HLL::calc_flux(const mmf::PState &zL, const mmf::PState &zR, const MixturePT &mix) {
    using namespace mmf;

    // Нормальные скорости слева и справа
    const double& u_L = zL.velocity.x();
    const double& u_R = zR.velocity.x();

    // Скорость звука слева и справа
    double c_L = mix.sound_speed_rP(zL.density, zL.pressure, zL.mass_frac,
                                    {.T0 = zL.T(), .rhos = &zL.densities});
    double c_R = mix.sound_speed_rP(zR.density, zR.pressure, zR.mass_frac,
                                    {.T0 = zR.T(), .rhos = &zR.densities});

    // Оценки скоростей расходящихся волн
    double S_L = min(u_L - c_L, u_R - c_R, 0.0);
    double S_R = max(u_L + c_L, u_R + c_R, 0.0);

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

mmf::WaveConfig2 HLL::wave_config(const MixturePT& mix,
                                  const mmf::PState& zL, const mmf::PState& zR) {
    // Нормальные скорости слева и справа
    const double &u_L = zL.velocity.x();
    const double &u_R = zR.velocity.x();

    // Скорость звука слева и справа
    double c_L = mix.sound_speed_rP(zL.density, zL.pressure, zL.mass_frac, {.T0=zL.T(), .rhos=&zL.rhos()});
    double c_R = mix.sound_speed_rP(zR.density, zR.pressure, zR.mass_frac, {.T0=zR.T(), .rhos=&zR.rhos()});

    // Оценки скоростей расходящихся волн
    double S_L = min(u_L - c_L, u_R - c_R, 0.0);
    double S_R = max(u_L + c_L, u_R + c_R, 0.0);

    mmf::QState Q_L(zL); // Консервативный вектор слева
    mmf::QState Q_R(zR); // Консервативный вектор справа

    mmf::Flux F_L(zL);   // Дифференциальный поток слева
    mmf::Flux F_R(zR);   // Дифференциальный поток справа

    mmf::Flux F = (S_R * F_L.arr() - S_L * F_R.arr() + S_L * S_R * (Q_R.arr() - Q_L.arr())) / (S_R - S_L);
    mmf::QState Q = (F_L.arr() - F_R.arr() + S_R * Q_R.arr() - S_L * Q_L.arr()) / (S_R - S_L);

    return mmf::WaveConfig2({.S_L = S_L, .S_R = S_R, .Qs = Q, .Fs = F});
}

mmf::WaveConfig2 HLL::wave_config(const MixturePT& mix,
                                  const mmf::QState& Q_L, const mmf::Flux& F_L,
                                  const mmf::QState& Q_R, const mmf::Flux& F_R) {
    Vector3d vL = Q_L.momentum / Q_L.density;
    Vector3d vR = Q_R.momentum / Q_R.density;

    double eL = Q_L.energy / Q_L.density - 0.5 * vL.squaredNorm();
    double eR = Q_R.energy / Q_R.density - 0.5 * vR.squaredNorm();

    Fractions beta_L = Q_L.mass_frac.arr() / Q_L.density;
    Fractions beta_R = Q_R.mass_frac.arr() / Q_R.density;

    // Скорость звука слева и справа
    double c_L = mix.sound_speed_re(Q_L.density, eL, beta_L);
    double c_R = mix.sound_speed_re(Q_R.density, eR, beta_R);

    // Оценки скоростей расходящихся волн
    double S_L = min(vL.x() - c_L, vR.x() - c_R, 0.0);
    double S_R = max(vL.x() + c_L, vR.x() + c_R, 0.0);

    // Поток через левую/правую волну
    mmf::Flux a_L = F_L.arr() - S_L * Q_L.arr();
    mmf::Flux a_R = F_R.arr() - S_R * Q_R.arr();

    mmf::Flux F = (S_R * a_L.arr() - S_L * a_R.arr()) / (S_R - S_L);
    mmf::QState Q = (a_L.arr() - a_R.arr()) / (S_R - S_L);

#if 0
    double err1 = ((F_R.arr() - F.arr()) - S_R * (Q_R.arr() - Q.arr())).cwiseAbs().maxCoeff();
    double err2 = ((F.arr() - F_L.arr()) - S_L * (Q.arr() - Q_L.arr())).cwiseAbs().maxCoeff();

    if (std::max({err1, err2}) > 1.0e-6) {
        throw std::runtime_error("Rankine-Hugoniot conditions failed");
    }
#endif

    return mmf::WaveConfig2({.S_L = S_L, .S_R = S_R, .Qs = Q, .Fs = F});
}

} // namespace zephyr::math