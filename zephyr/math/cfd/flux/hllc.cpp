#include <iostream>

#include <zephyr/geom/vector.h>
#include <zephyr/math/cfd/flux/hllc.h>

namespace zephyr::math {

smf::Flux HLLC2::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux HLLC2::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    double rho1 = zL.density, rho2 = zR.density;
    double u1 = zL.velocity.x(), u2 = zR.velocity.x();
    double p1 = zL.pressure, p2 = zR.pressure;

    double c1 = eos.sound_speed_rp(rho1, p1); // скорость звука слева std::sqrt(gamma * pressure / density);
    double c2 = eos.sound_speed_rp(rho2, p2); // скорость звука справа

    double s1 = std::min(u1 - c1, u2 - c2); // SL
    double s2 = std::max(u1 + c1, u2 + c2); // SR

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа
    if (0 < s1)
        return fL;
    if (s2 < 0)
        return fR;

    double S = (p2 - p1 + rho1 * u1 * (s1 - u1) - rho2 * u2 * (s2 - u2)) / (rho1 * (s1 - u1) - rho2 * (s2 - u2)); // S*

    double rhoL = rho1 * (s1 - u1) / (s1 - S); // rho*_L
    double rhoR = rho2 * (s2 - u2) / (s2 - S); // rho*_R

    double P = p1 + rho1 * (s1 - u1) * (S - u1); // P*

    Flux F_hllc;
    if (s1 <= 0 && 0 <= S) {
        QState qL(zL); // U_L
        double EL = ((s1 - u1) * qL.energy + S * P - u1 * p1) / (s1 - S);
        QState QL(rhoL, rhoL * Vector3d(S, zL.velocity.y(), zL.velocity.z()), EL); // U*_L

        F_hllc = fL.vec() + s1 * (QL.vec() - qL.vec());
    } else if (S <= 0 && 0 <= s2) {
        QState qR(zR); // U_R
        double ER = ((s2 - u2) * qR.energy + S * P - u2 * p2) / (s2 - S);
        QState QR(rhoR, rhoR * Vector3d(S, zR.velocity.y(), zR.velocity.z()), ER); // U*_R

        F_hllc = fR.vec() + s2 * (QR.vec() - qR.vec());
    } else {
        std::cerr << "zL: " << zL << "\n"; // сюда попадает только если где-то none или другие кривые значения
        std::cerr << "Zr: " << zR << "\n";
        std::cerr << "Sound speed left: " << c1 << " , Sound speed right: " << c2 << "\n";
        std::cerr << "SL: " << s1 << ", SR: " << s2 << ", S*: " << S << "\n";
        throw std::runtime_error("HLLC2::calc_flux Error, strange case in switch");
    }

    return F_hllc;
}

mmf::Flux HLLC2::mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const {
    using namespace mmf;

    double rho1 = zL.density, rho2 = zR.density;
    double u1 = zL.velocity.x(), u2 = zR.velocity.x();
    double p1 = zL.pressure, p2 = zR.pressure;

    double c1 = mixture.sound_speed_rp(rho1, p1, zL.mass_frac); // скорость звука слева
    double c2 = mixture.sound_speed_rp(rho2, p2, zR.mass_frac); // скорость звука справа

    double s1 = std::min(u1 - c1, u2 - c2); // SL
    double s2 = std::max(u1 + c1, u2 + c2); // SR

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа
    if (0 < s1)
        return fL;
    if (s2 < 0)
        return fR;

    double S = (p2 - p1 + rho1 * u1 * (s1 - u1) - rho2 * u2 * (s2 - u2)) / (rho1 * (s1 - u1) - rho2 * (s2 - u2)); // S*

    double rhoL = rho1 * (s1 - u1) / (s1 - S); // rho*_L
    double rhoR = rho2 * (s2 - u2) / (s2 - S); // rho*_R

    double P = p1 + rho1 * (s1 - u1) * (S - u1); // P*

    Flux F_hllc;
    if (s1 <= 0 && 0 <= S) {
        QState qL(zL); // U_L
        double EL = ((s1 - u1) * qL.energy + S * P - u1 * p1) / (s1 - S);
        // TODO: здесь надо как-то рассчитать mass_frac
        QState QL(rhoL, rhoL * Vector3d(S, zL.velocity.y(), zL.velocity.z()), EL, FractionsFlux(zL.mass_frac)); // U*_L

        F_hllc = fL.vec() + s1 * (QL.vec() - qL.vec());
    } else if (S <= 0 && 0 <= s2) {
        QState qR(zR); // U_R
        double ER = ((s2 - u2) * qR.energy + S * P - u2 * p2) / (s2 - S);
        // TODO: здесь надо как-то рассчитать mass_frac
        QState QR(rhoR, rhoR * Vector3d(S, zR.velocity.y(), zR.velocity.z()), ER, FractionsFlux(zR.mass_frac)); // U*_R

        F_hllc = fR.vec() + s2 * (QR.vec() - qR.vec());
    } else {
        std::cerr << "zL: " << zL << "\n"; // сюда попадает только если где-то none или другие кривые значения
        std::cerr << "Zr: " << zR << "\n";
        std::cerr << "Sound speed left: " << c1 << " , Sound speed right: " << c2 << "\n";
        std::cerr << "SL: " << s1 << ", SR: " << s2 << ", S*: " << S << "\n";
        throw std::runtime_error("HLLC2::calc_flux Error, strange case in switch");
    }

    return F_hllc;
}

smf::Flux HLLC::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux HLLC::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    double rho1 = zL.density, rho2 = zR.density;
    double u1 = zL.velocity.x(), u2 = zR.velocity.x();
    double v1 = zL.velocity.y(), v2 = zR.velocity.y();
    double w1 = zL.velocity.z(), w2 = zR.velocity.z();
    double p1 = zL.pressure, p2 = zR.pressure;

    QState qL(zL); // U_L
    QState qR(zR); // U_R
    double E1 = qL.energy, E2 = qR.energy;

    double c1 = eos.sound_speed_rp(rho1, p1); // скорость звука слева std::sqrt(gamma * pressure / density);
    double c2 = eos.sound_speed_rp(rho2, p2); // скорость звука справа

    double s1 = std::min(u1 - c1, u2 - c2); // SL
    double s2 = std::max(u1 + c1, u2 + c2); // SR
    double S = (p2 - p1 + rho1 * u1 * (s1 - u1) - rho2 * u2 * (s2 - u2)) / (rho1 * (s1 - u1) - rho2 * (s2 - u2)); // S*

    double rhoL = rho1 * (s1 - u1) / (s1 - S); // rho*_L
    double rhoR = rho2 * (s2 - u2) / (s2 - S); // rho*_R

    double P = p1 + rho1 * (s1 - u1) * (S - u1); // P*

    double EL = ((s1 - u1) * E1 + S * P - u1 * p1) / (s1 - S); // сам вывел
    double ER = ((s2 - u2) * E2 + S * P - u2 * p2) / (s2 - S);
    
    // double EL = rhoL * (E1 / rho1 + (S - u1) * (S + p1 / (rho1 * (s1 - u1)))); // из Тора
    // double ER = rhoR * (E2 / rho2 + (S - u2) * (S + p2 / (rho2 * (s2 - u2))));

    QState QL(rhoL, rhoL * Vector3d(S, v1, w1), EL); // U*_L
    QState QR(rhoR, rhoR * Vector3d(S, v2, w2), ER); // U*_R

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа

    Flux F_hllc;
    if (0 < s1)
        F_hllc = fL;
    else if (s1 <= 0 && 0 <= S)
        F_hllc = fL.vec() + s1 * (QL.vec() - qL.vec());
    else if (S <= 0 && 0 <= s2)
        F_hllc = fR.vec() + s2 * (QR.vec() - qR.vec());
    else if (s2 < 0)
        F_hllc = fR;
    else {
        std::cerr << "zL: " << zL << "\n"; // сюда попадает только если где-то none или другие кривые значения
        std::cerr << "Zr: " << zR << "\n";
        std::cerr << "Sound speed left: " << c1 << " , Sound speed right: " << c2 << "\n";
        std::cerr << "SL: " << s1 << ", SR: " << s2 << ", S*: " << S << "\n";
        throw std::runtime_error("HLLC::calc_flux Error, strange case in switch");
    }

    /*
    if (std::isnan(F_hllc.mass) || std::isnan(F_hllc.momentum.x()) || std::isnan(F_hllc.momentum.y()) ||
        std::isnan(F_hllc.momentum.z()) || std::isnan(F_hllc.energy)) {
        std::cerr << "zL: " << zL << "\n";
        std::cerr << "Zr: " << zR << "\n";
        std::cerr << "Sound speed left: " << c1 << " , Sound speed right: " << c2 << "\n";
        std::cerr << "SL: " << s1 << ", SR: " << s2 << ", S*: " << S << "\n";
        std::cerr << "F_HLLC: " << F_hllc << "\n";
        throw std::runtime_error("HLLC::calc_flux Error, F_HLLC has the bad value");
    }
     */

    return F_hllc;
}

mmf::Flux HLLC::mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const {
    using namespace mmf;

    double rho1 = zL.density, rho2 = zR.density;
    double u1 = zL.velocity.x(), u2 = zR.velocity.x();
    double v1 = zL.velocity.y(), v2 = zR.velocity.y();
    double w1 = zL.velocity.z(), w2 = zR.velocity.z();
    double p1 = zL.pressure, p2 = zR.pressure;

    QState qL(zL); // U_L
    QState qR(zR); // U_R
    double E1 = qL.energy, E2 = qR.energy;

    double c1 = mixture.sound_speed_rp(rho1, p1, zL.mass_frac); // скорость звука слева
    double c2 = mixture.sound_speed_rp(rho2, p2, zR.mass_frac); // скорость звука справа

    double s1 = std::min(u1 - c1, u2 - c2); // SL
    double s2 = std::max(u1 + c1, u2 + c2); // SR
    double S = (p2 - p1 + rho1 * u1 * (s1 - u1) - rho2 * u2 * (s2 - u2)) / (rho1 * (s1 - u1) - rho2 * (s2 - u2)); // S*

    double rhoL = rho1 * (s1 - u1) / (s1 - S); // rho*_L
    double rhoR = rho2 * (s2 - u2) / (s2 - S); // rho*_R

    double P = p1 + rho1 * (s1 - u1) * (S - u1); // P*

    double EL = ((s1 - u1) * E1 + S * P - u1 * p1) / (s1 - S);
    double ER = ((s2 - u2) * E2 + S * P - u2 * p2) / (s2 - S);

    // TODO: здесь надо как-то рассчитать mass_frac
    QState QL(rhoL, rhoL * Vector3d(S, v1, w1), EL, FractionsFlux(zL.mass_frac)); // U*_L
    QState QR(rhoR, rhoR * Vector3d(S, v2, w2), ER, FractionsFlux(zR.mass_frac)); // U*_R

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа

    Flux F_hllc;
    if (0 < s1)
        F_hllc = fL;
    else if (s1 <= 0 && 0 <= S)
        F_hllc = fL.vec() + s1 * (QL.vec() - qL.vec());
    else if (S <= 0 && 0 <= s2)
        F_hllc = fR.vec() + s2 * (QR.vec() - qR.vec());
    else if (s2 < 0)
        F_hllc = fR;
    else {
        std::cerr << "zL: " << zL << "\n"; // сюда попадает только если где-то none или другие кривые значения
        std::cerr << "Zr: " << zR << "\n";
        std::cerr << "Sound speed left: " << c1 << " , Sound speed right: " << c2 << "\n";
        std::cerr << "SL: " << s1 << ", SR: " << s2 << ", S*: " << S << "\n";
        throw std::runtime_error("HLLC::calc_flux Error, strange case in switch");
    }

    return F_hllc;
}

smf::Flux HLLC_central::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux HLLC_central::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    double rho1 = zL.density, rho2 = zR.density;
    double u1 = zL.velocity.x(), u2 = zR.velocity.x();
    double v1 = zL.velocity.y(), v2 = zR.velocity.y();
    double w1 = zL.velocity.z(), w2 = zR.velocity.z();
    double p1 = zL.pressure, p2 = zR.pressure;

    QState qL(zL); // U_L
    QState qR(zR); // U_R
    double E1 = qL.energy, E2 = qR.energy;

    double c1 = eos.sound_speed_rp(rho1, p1); // скорость звука слева std::sqrt(gamma * pressure / density);
    double c2 = eos.sound_speed_rp(rho2, p2); // скорость звука справа

    // uˆ and cˆ are determined from the Roe average
    double u = (u1 * sqrt(rho1) + u2 * sqrt(rho2)) / (sqrt(rho1) + sqrt(rho2));
    double c = sqrt((std::pow(c1, 2) * sqrt(rho1) + std::pow(c2, 2) * sqrt(rho2)) / (sqrt(rho1) + sqrt(rho2)) + 0.5 * sqrt(rho1 * rho2) * std::pow(u1 - u2, 2) / std::pow(sqrt(rho1) + sqrt(rho2), 2));

    double s1 = std::min(u1 - c1, u - c); // SL
    double s2 = std::max(u1 + c1, u + c); // SR
    double S = (p2 - p1 + rho1 * u1 * (s1 - u1) - rho2 * u2 * (s2 - u2)) / (rho1 * (s1 - u1) - rho2 * (s2 - u2)); // S*

    double k1 = (s1 - u1) / (s1 - S);
    double k2 = (s2 - u2) / (s2 - S);

    double mass1 = k1 * rho1;
    double mass2 = k2 * rho2;

    Vector3d momentum1 = k1 * rho1 * Vector3d(S, v1, w1);
    Vector3d momentum2 = k2 * rho2 * Vector3d(S, v2, w2);

    double energy1 = k1 * (E1 + (S - u1) * (rho1 * S + p1 / (s1 - S)));
    double energy2 = k2 * (E2 + (S - u2) * (rho2 * S + p2 / (s2 - S)));

    QState QL(mass1, momentum1, energy1);
    QState QR(mass2, momentum2, energy2);

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа

    if (s1 >= 0)
        return fL;
    if (s2 <= 0)
        return fR;
    return 0.5 * (fL.vec() + fR.vec()) + 0.5 * (s1 * (QL.vec() - qL.vec()) + std::abs(S) * (QL.vec() - QR.vec()) + s2 * (QR.vec() - qR.vec()));
}

mmf::Flux HLLC_central::mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const {}

smf::Flux HLLC_LM::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux HLLC_LM::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    double rho1 = zL.density, rho2 = zR.density;
    double u1 = zL.velocity.x(), u2 = zR.velocity.x();
    double v1 = zL.velocity.y(), v2 = zR.velocity.y();
    double w1 = zL.velocity.z(), w2 = zR.velocity.z();
    double p1 = zL.pressure, p2 = zR.pressure;
    
    QState qL(zL); // U_L
    QState qR(zR); // U_R
    double E1 = qL.energy, E2 = qR.energy;

    double c1 = eos.sound_speed_rp(rho1, p1); // скорость звука слева std::sqrt(gamma * pressure / density);
    double c2 = eos.sound_speed_rp(rho2, p2); // скорость звука справа

    double Malimit = 0.1;
    double Malocal = std::max(std::abs(u1/c1), std::abs(u2/c2));
    double phi = std::sin(std::min(1.0, Malocal / Malimit) * M_PI_2);

    // uˆ and cˆ are determined from the Roe average
    double u = (u1 * sqrt(rho1) + u2 * sqrt(rho2)) / (sqrt(rho1) + sqrt(rho2));
    double c = sqrt((std::pow(c1, 2) * sqrt(rho1) + std::pow(c2, 2) * sqrt(rho2)) / (sqrt(rho1) + sqrt(rho2)) + 0.5 * sqrt(rho1 * rho2) * std::pow(u1 - u2, 2) / std::pow(sqrt(rho1) + sqrt(rho2), 2));

    double s1 = std::min(u1 - c1, u - c); // SL
    double s2 = std::max(u1 + c1, u + c); // SR
    double S = (p2 - p1 + rho1 * u1 * (s1 - u1) - rho2 * u2 * (s2 - u2)) / (rho1 * (s1 - u1) - rho2 * (s2 - u2)); // S*

    double k1 = (s1 - u1) / (s1 - S);
    double k2 = (s2 - u2) / (s2 - S);

    double mass1 = k1 * rho1;
    double mass2 = k2 * rho2;

    Vector3d momentum1 = k1 * rho1 * Vector3d(S, v1, w1);
    Vector3d momentum2 = k2 * rho2 * Vector3d(S, v2, w2);

    double energy1 = k1 * (E1 + (S - u1) * (rho1 * S + p1 / (s1 - S)));
    double energy2 = k2 * (E2 + (S - u2) * (rho2 * S + p2 / (s2 - S)));

    QState QL(mass1, momentum1, energy1);
    QState QR(mass2, momentum2, energy2);

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа

    if (s1 >= 0)
        return fL;
    if (s2 <= 0)
        return fR;
    return 0.5 * (fL.vec() + fR.vec()) + 0.5 * (phi * s1 * (QL.vec() - qL.vec()) + std::abs(S) * (QL.vec() - QR.vec()) + phi * s2 * (QR.vec() - qR.vec()));
}

mmf::Flux HLLC_LM::mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const {}

}