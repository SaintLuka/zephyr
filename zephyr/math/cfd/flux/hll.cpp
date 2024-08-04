#include <iostream>

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
    double c_L = eos.sound_speed_rp(zL.density, zL.pressure);
    double c_R = eos.sound_speed_rp(zR.density, zR.pressure);

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

mmf::Flux HLL::flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const {
    return calc_flux(zL, zR, mixture);
}

mmf::Flux HLL::calc_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) {
    using namespace mmf;

    // Нормальные скорости слева и справа
    const double& u_L = zL.velocity.x();
    const double& u_R = zR.velocity.x();

    // Скорость звука слева и справа
    double c_L = mixture.sound_speed_rp(zL.density, zL.pressure,
                                        zL.mass_frac, {.T0 = zL.temperature});
    double c_R = mixture.sound_speed_rp(zR.density, zR.pressure,
                                        zR.mass_frac, {.T0 = zR.temperature});

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