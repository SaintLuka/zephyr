#include <iostream>

#include <zephyr/geom/vector.h>
#include <zephyr/math/cfd/flux/rusanov.h>

namespace zephyr::math {

smf::Flux Rusanov::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux Rusanov::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
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

    // Максимальное собственное значение
    double S = std::max(std::abs(S_L), S_R);

    QState Q_L(zL); // Консервативный вектор слева
    QState Q_R(zR); // Консервативный вектор справа

    Flux F_L(zL);   // Дифференциальный поток слева
    Flux F_R(zR);   // Дифференциальный поток справа

    Flux F = 0.5 * (F_L.vec() + F_R.vec()) + 0.5 * S * (Q_L.vec() - Q_R.vec());

    return F;
}

mmf::Flux Rusanov::flux(const mmf::PState &zL, const mmf::PState &zR, const phys::MixturePT &mixture) const {
    return calc_flux(zL, zR, mixture);
}

mmf::Flux Rusanov::calc_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::MixturePT &mixture) {
    using namespace mmf;

    // Нормальные скорости слева и справа
    const double &u_L = zL.velocity.x();
    const double &u_R = zR.velocity.x();

    // Скорость звука слева и справа
    double c_L = mixture.sound_speed_rP(zL.density, zL.pressure, zL.mass_frac,
                                        {.T0 = zL.T(), .rhos = &zL.densities});
    double c_R = mixture.sound_speed_rP(zR.density, zR.pressure, zR.mass_frac,
                                        {.T0 = zR.T(), .rhos = &zR.densities});

    if (std::isnan(c_L) || std::isnan(c_R)) {
        std::cerr << "Rusanov::calc_mm_flux error #1\n";
        std::cerr << "Can't find sound speed for state: " << zL << "\n";
        std::cerr << "Can't find sound speed for state: " << zR << "\n";
        throw std::runtime_error("Rusanov::calc_mm_flux error #1");
    }

    // Оценки скоростей расходящихся волн
    double S_L = std::min({u_L - c_L, u_R - c_R, 0.0});
    double S_R = std::max({u_L + c_L, u_R + c_R, 0.0});

    // Максимальное собственное значение
    double S = std::max(std::abs(S_L), S_R);

    QState Q_L(zL); // Консервативный вектор слева
    QState Q_R(zR); // Консервативный вектор справа

    Flux F_L(zL);   // Дифференциальный поток слева
    Flux F_R(zR);   // Дифференциальный поток справа

    Flux F = 0.5 * (F_L.vec() + F_R.vec()) + 0.5 * S * (Q_L.vec() - Q_R.vec());

    if (F.is_bad()) {
        std::cerr << "Rusanov::calc_mm_flux error #2\n";
        std::cerr << "  zL: " << zL << "; c_L: " << c_L << "\n";
        std::cerr << "  zR: " << zR << "; c_R: " << c_R << "\n";
        std::cerr << "  Flux: " << F << "\n";
        throw std::runtime_error("Rusanov::calc_mm_flux error #2");
    }

    return F;
}

} // namespace zephyr::math