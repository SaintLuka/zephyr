#include <iostream>
#include <iomanip>

#include <zephyr/geom/vector.h>
#include <zephyr/math/cfd/flux/hllc.h>

namespace zephyr::math {

inline double sqr(double x) {
    return x * x;
}

smf::Flux HLLC::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux HLLC::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    const double &rho_L = zL.density;
    const double &rho_R = zR.density;
    const double &u_L = zL.velocity.x();
    const double &u_R = zR.velocity.x();
    const double &P_L = zL.pressure;
    const double &P_R = zR.pressure;

    // Скорость звука слева и справа
    double c_L = eos.sound_speed_rp(rho_L, P_L);
    double c_R = eos.sound_speed_rp(rho_R, P_R);

    // Оценки скоростей расходящихся волн
    double S_L = std::min({u_L - c_L, u_R - c_R, 0.0});
    double S_R = std::max({u_L + c_L, u_R + c_R, 0.0});

    // Скорость контактного разрыва
    double S_C = (P_R - P_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)) /
                 (rho_L * (S_L - u_L) - rho_R * (S_R - u_R));

    // Давление на контактном разрыве
    double P_C = P_L + rho_L * (S_L - u_L) * (S_C - u_L);

    Flux F;
    if (S_C >= 0.0) {
        // Плотность слева от контактного разрыва
        double rho_sL = rho_L * (S_L - u_L) / (S_L - S_C);

        // Внутренняя энергия слева от контактного разрыва
        double e_sL = zL.energy + 0.5 * sqr(S_C - u_L) + P_L * (S_C - u_L) / (rho_L * (S_L - u_L));

        // Вектор состояния слева от контактного разрыва
        PState z_sL(rho_sL, {S_C, zL.v(), zL.w()}, P_C, e_sL);

        // Консервативный вектор слева от контактного разрыва
        QState Q_sL(z_sL);

        QState Q_L(zL);  // Консервативный вектор слева
        Flux F_L(zL);  // Дифференциальный поток слева

        F = F_L.arr() + S_L * (Q_sL.arr() - Q_L.arr());
    } else {
        // Плотность слева от контактного разрыва
        double rho_sR = rho_R * (S_R - u_R) / (S_R - S_C);

        // Внутренняя энергия слева от контактного разрыва
        double e_sR = zR.energy + 0.5 * sqr(S_C - u_R) + P_R * (S_C - u_R) / (rho_R * (S_R - u_R));

        // Вектор состояния слева от контактного разрыва
        PState z_sR(rho_sR, {S_C, zR.v(), zR.w()}, P_C, e_sR);

        // Консервативный вектор слева от контактного разрыва
        QState Q_sR(z_sR);

        QState Q_R(zR);  // Консервативный вектор слева
        Flux F_R(zR);    // Дифференциальный поток справа

        F = F_R.arr() + S_R * (Q_sR.arr() - Q_R.arr());
    }

    if (F.arr().hasNaN()) {
        std::cerr << "HLLC::calc_flux error\n";
        std::cerr << "  z_L: " << zL << "\n";
        std::cerr << "  z_R: " << zR << "\n";
        std::cerr << "  c_L: " << c_L << "; c_R: " << c_R << "\n";
        std::cerr << "  S_L: " << S_L << "; S_R: " << S_R << "\n";
        std::cerr << "  F_HLL: " << F << "\n";
        throw std::runtime_error("HLLC::calc_flux error: bad value");
    }

    return F;
}

mmf::Flux HLLC::flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const {
    return calc_flux(zL, zR, mixture);
}

mmf::Flux HLLC::calc_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) {
    using namespace mmf;

    const double &rho_L = zL.density;
    const double &rho_R = zR.density;
    const double &u_L = zL.velocity.x();
    const double &u_R = zR.velocity.x();
    const double &P_L = zL.pressure;
    const double &P_R = zR.pressure;

    // Скорость звука слева и справа
    double c_L = mixture.sound_speed_rp(zL.density, zL.pressure,
                                        zL.mass_frac, {.T0 = zL.temperature});
    double c_R = mixture.sound_speed_rp(zR.density, zR.pressure,
                                        zR.mass_frac, {.T0 = zR.temperature});

    // Оценки скоростей расходящихся волн
    double S_L = std::min({u_L - c_L, u_R - c_R, 0.0});
    double S_R = std::max({u_L + c_L, u_R + c_R, 0.0});

    // Скорость контактного разрыва
    double S_C = (P_R - P_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)) /
                 (rho_L * (S_L - u_L) - rho_R * (S_R - u_R));

    // Давление на контактном разрыве
    double P_C = P_L + rho_L * (S_L - u_L) * (S_C - u_L);

    Flux F;
    if (S_C >= 0.0) {
        // Плотность слева от контактного разрыва
        double rho_sL = rho_L * (S_L - u_L) / (S_L - S_C);

        // Внутренняя энергия слева от контактного разрыва
        double e_sL = zL.energy + 0.5 * sqr(S_C - u_L) + P_L * (S_C - u_L) / (rho_L * (S_L - u_L));

        // Вектор состояния слева от контактного разрыва
        PState z_sL(rho_sL, {S_C, zL.v(), zL.w()}, P_C, e_sL, zL.beta(), NAN, Fractions::NaN());

        // Консервативный вектор слева от контактного разрыва
        QState Q_sL(z_sL);

        QState Q_L(zL);  // Консервативный вектор слева
        Flux F_L(zL);    // Дифференциальный поток слева

        F = F_L.arr() + S_L * (Q_sL.arr() - Q_L.arr());
    } else {
        // Плотность слева от контактного разрыва
        double rho_sR = rho_R * (S_R - u_R) / (S_R - S_C);

        // Внутренняя энергия слева от контактного разрыва
        double e_sR = zR.energy + 0.5 * sqr(S_C - u_R) + P_R * (S_C - u_R) / (rho_R * (S_R - u_R));

        // Вектор состояния слева от контактного разрыва
        PState z_sR(rho_sR, {S_C, zR.v(), zR.w()}, P_C, e_sR, zR.beta(), NAN, Fractions::NaN());

        // Консервативный вектор слева от контактного разрыва
        QState Q_sR(z_sR);

        QState Q_R(zR);  // Консервативный вектор слева
        Flux F_R(zR);    // Дифференциальный поток справа

        F = F_R.arr() + S_R * (Q_sR.arr() - Q_R.arr());
    }

    if (F.arr().hasNaN()) {
        std::cerr << "HLLC::calc_flux error\n";
        std::cerr << "  z_L: " << zL << "\n";
        std::cerr << "  z_R: " << zR << "\n";
        std::cerr << "  c_L: " << c_L << "; c_R: " << c_R << "\n";
        std::cerr << "  S_L: " << S_L << "; S_R: " << S_R << "\n";
        std::cerr << "  F_HLL: " << F << "\n";
        throw std::runtime_error("HLLC::calc_flux error: bad value");
    }

    return F;

}

smf::Flux HLLC_C::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux HLLC_C::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
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
    return 0.5 * (fL.arr() + fR.arr()) + 0.5 * (s1 * (QL.arr() - qL.arr()) + std::abs(S) * (QL.arr() - QR.arr()) + s2 * (QR.arr() - qR.arr()));
}

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
    return 0.5 * (fL.arr() + fR.arr()) + 0.5 * (phi * s1 * (QL.arr() - qL.arr()) + std::abs(S) * (QL.arr() - QR.arr()) + phi * s2 * (QR.arr() - qR.arr()));
}

}