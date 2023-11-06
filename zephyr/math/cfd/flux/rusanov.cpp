#include <iostream>

#include <zephyr/geom/vector.h>
#include <zephyr/math/cfd/flux/rusanov.h>

namespace zephyr::math {

smf::Flux Rusanov::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux Rusanov::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    // Состояние на грани
    PState zE = 0.5 * (zL.vec() + zR.vec());

    // Согласуем состояние на грани
    zE.energy = eos.energy_rp(zE.density, zE.pressure);

    double rho = zE.density;
    double u = zE.velocity.x();
    double p = zE.pressure;

    double c = eos.sound_speed_rp(rho, p);

    // Максимальное собственное значение
    double max_lambda = std::abs(u) + c;

    QState qL(zL); // Консервативный вектор слева
    QState qR(zR); // Консервативный вектор справа

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа

    Flux res = 0.5 * (fL.vec() + fR.vec()) + 0.5 * max_lambda * (qL.vec() - qR.vec());

    return res;
}

smf::Flux Rusanov2::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux Rusanov2::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    double rho1 = zL.density, rho2 = zR.density;
    double u1 = zL.velocity.x(), u2 = zR.velocity.x();
    double p1 = zL.pressure, p2 = zR.pressure;

    double c1 = eos.sound_speed_rp(rho1, p1), c2 = eos.sound_speed_rp(rho2, p2);

    // Максимальное собственное значение
    double max_lambda = std::max(std::abs(u1) + c1, std::abs(u2) + c2);

    QState qL(zL); // Консервативный вектор слева
    QState qR(zR); // Консервативный вектор справа

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа

    Flux res = 0.5 * (fL.vec() + fR.vec()) + 0.5 * max_lambda * (qL.vec() - qR.vec());

    return res;
}

mmf::Flux Rusanov2::mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const {
    using namespace mmf;

    double rho1 = zL.density, rho2 = zR.density;
    double u1 = zL.velocity.x(), u2 = zR.velocity.x();
    double p1 = zL.pressure, p2 = zR.pressure;

    double c1 = mixture.sound_speed_rp(rho1, p1, zL.mass_frac);
    double c2 = mixture.sound_speed_rp(rho2, p2, zR.mass_frac);

    // Максимальное собственное значение
    double max_lambda = std::max(std::abs(u1) + c1, std::abs(u2) + c2);

    QState qL(zL); // Консервативный вектор слева
    QState qR(zR); // Консервативный вектор справа

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа

    Flux res = 0.5 * (fL.vec() + fR.vec()) + 0.5 * max_lambda * (qL.vec() - qR.vec());

    return res;
}


} // namespace zephyr