#include <iostream>

#include <zephyr/geom/vector.h>
#include <zephyr/math/cfd/flux/hll.h>

namespace zephyr::math {

smf::Flux HLL::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux HLL::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    double rho1 = zL.density, rho2 = zR.density;
    double u1 = zL.velocity.x(), u2 = zR.velocity.x();
    double p1 = zL.pressure, p2 = zR.pressure;

    QState qL(zL); // U_L
    QState qR(zR); // U_R

    double c1 = eos.sound_speed_rp(rho1, p1); // скорость звука слева std::sqrt(gamma * pressure / density);
    double c2 = eos.sound_speed_rp(rho2, p2); // скорость звука справа

    double s1 = std::min(u1 - c1, u2 - c2); // SL
    double s2 = std::max(u1 + c1, u2 + c2); // SR

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа

    Flux F_hll;
    if (0 < s1)
        F_hll = fL;
    else if (s1 <= 0 && 0 <= s2)
        F_hll = (s2 * fL.vec() - s1 * fR.vec() + s1 * s2 * (qR.vec() - qL.vec())) / (s2 - s1);
    else if (s2 < 0)
        F_hll = fR;
    else {
        std::cerr << "zL: " << zL << "\n"; // сюда попадает только если где-то none или другие кривые значения
        std::cerr << "Zr: " << zR << "\n";
        std::cerr << "Sound speed left: " << c1 << " , Sound speed right: " << c2 << "\n";
        std::cerr << "SL: " << s1;
        throw std::runtime_error("HLL::calc_flux Error, strange case in switch");
    }

    if (std::isnan(F_hll.mass) || std::isnan(F_hll.momentum.x()) || std::isnan(F_hll.momentum.y()) ||
        std::isnan(F_hll.momentum.z()) || std::isnan(F_hll.energy)) {
        std::cerr << "zL: " << zL << "\n";
        std::cerr << "Zr: " << zR << "\n";
        std::cerr << "Sound speed left: " << c1 << " , Sound speed right: " << c2 << "\n";
        std::cerr << "SL: " << s1 << ", SR: " << s2;
        std::cerr << "F_HLL: " << F_hll << "\n";
        throw std::runtime_error("HLL::calc_flux Error, F_HLL has the bad value");
    }

    return F_hll;
}

mmf::Flux HLL::mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const {
    return HLL::calc_mm_flux(zL, zR, mixture);
}

mmf::Flux HLL::calc_mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) {
    using namespace mmf;

    QState qL(zL); // U_L
    QState qR(zR); // U_R

    double c1 = mixture.sound_speed_rp(zL.density, zL.pressure, zL.mass_frac, {.T0 = zL.temperature}); // скорость звука слева
    double c2 = mixture.sound_speed_rp(zR.density, zR.pressure, zR.mass_frac, {.T0 = zR.temperature}); // скорость звука справа

    double u1 = zL.velocity.x(), u2 = zR.velocity.x();
    double s1 = std::min({u1 - c1, u2 - c2, 0.0}); // SL
    double s2 = std::max({u1 + c1, u2 + c2, 0.0}); // SR

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа

    Flux F_hll;
    if (0 < s1)
        F_hll = fL;
    else if (s1 <= 0 && 0 <= s2)
        F_hll = (s2 * fL.vec() - s1 * fR.vec() + s1 * s2 * (qR.vec() - qL.vec())) / (s2 - s1);
    else if (s2 < 0)
        F_hll = fR;
    else {
        std::cerr << "zL: " << zL << "\n"; // сюда попадает только если где-то none или другие кривые значения
        std::cerr << "Zr: " << zR << "\n";
        std::cerr << "Sound speed left: " << c1 << " , Sound speed right: " << c2 << "\n";
        std::cerr << "SL: " << s1;
        exit(1);
        throw std::runtime_error("HLL::calc_mm_flux Error, strange case in switch");
    }

    if (F_hll.is_bad()) {
        std::cerr << "zL: " << zL << "\n";
        std::cerr << "Zr: " << zR << "\n";
        std::cerr << "Sound speed left: " << c1 << " , Sound speed right: " << c2 << "\n";
        std::cerr << "SL: " << s1 << ", SR: " << s2;
        std::cerr << "F_HLL: " << F_hll << "\n";
        exit(1);
        throw std::runtime_error("HLL::calc_mm_flux Error, F_HLL has the bad value");
    }

    return F_hll;
}

}