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
    double c_L = eos.sound_speed_rP(rho_L, P_L);
    double c_R = eos.sound_speed_rP(rho_R, P_R);

    // Оценки скоростей расходящихся волн
    double S_L = std::min(u_L - c_L, u_R - c_R);
    double S_R = std::max(u_L + c_L, u_R + c_R);

    // Сверхзвуковое течение
    if (S_L >= 0.0) { return Flux(zL); }
    if (S_R <= 0.0) { return Flux(zR); }

    // Перенос массы через левую/правую волну
    double a_L = rho_L * (S_L - u_L);
    double a_R = rho_R * (S_R - u_R);

    // Скорость контактного разрыва
    double S_C = (P_L - P_R + a_R * u_R - a_L * u_L) / (a_R - a_L);

    // Плотность слева/справа от контактного разрыва
    double rho_sL = a_L / (S_L - S_C);
    double rho_sR = a_R / (S_R - S_C);

    // Удельная полная энергия слева/справа от контактного разрыва
    double E_sL = zL.E() + (S_C - u_L) * (S_C + P_L / a_L);
    double E_sR = zR.E() + (S_C - u_R) * (S_C + P_R / a_R);

    // Консервативный вектор слева/справа от контактного разрыва
    QState Q_sL(rho_sL, {rho_sL * S_C, rho_sL * zL.v(), rho_sL * zL.w()}, rho_sL * E_sL);
    QState Q_sR(rho_sR, {rho_sR * S_C, rho_sR * zR.v(), rho_sR * zR.w()}, rho_sR * E_sR);

    // Консервативный вектор слева/справа
    QState Q_L(zL);
    QState Q_R(zR);

    // Дифференциальный поток слева/справа
    Flux F_L(zL);
    Flux F_R(zR);

    // Симметричная формула HLL(C)
    Flux F = 0.5 * (F_L.arr() + F_R.arr() +
            S_L * (Q_sL.arr() - Q_L.arr()) +
            std::abs(S_C) * (Q_sL.arr() - Q_sR.arr()) +
            S_R * (Q_sR.arr() - Q_R.arr())
    );

    if (F.arr().hasNaN()) {
        std::cerr << "HLLC::calc_flux error\n";
        std::cerr << "  z_L: " << zL << "\n";
        std::cerr << "  z_R: " << zR << "\n";
        std::cerr << "  c_L: " << c_L << "; c_R: " << c_R << "\n";
        std::cerr << "  S_L: " << S_L << "; S_R: " << S_R << "\n";
        std::cerr << "  F_HLLC: " << F << "\n";
        throw std::runtime_error("HLLC::calc_flux error: bad value");
    }

    return F;
}

HLLC::WaveConfig HLLC::wave_config(
        const phys::Eos& eosL, const smf::PState& zL,
        const phys::Eos& eosR, const smf::PState& zR) {
    const double &rho_L = zL.density;
    const double &rho_R = zR.density;
    const double &u_L = zL.velocity.x();
    const double &u_R = zR.velocity.x();
    const double &P_L = zL.pressure;
    const double &P_R = zR.pressure;

    // Скорость звука слева и справа
    double c_L = eosL.sound_speed_rP(rho_L, P_L);
    double c_R = eosR.sound_speed_rP(rho_R, P_R);

    // Оценки скоростей расходящихся волн
    double S_L = std::min({u_L - c_L, u_R - c_R, 0.0});
    double S_R = std::max({u_L + c_L, u_R + c_R, 0.0});

    // Перенос массы через левую/правую волну
    double a_L = rho_L * (S_L - u_L);
    double a_R = rho_R * (S_R - u_R);

    // Скорость контактного разрыва
    double S_C = (P_L - P_R + a_R * u_R - a_L * u_L) / (a_R - a_L);

    // Плотность слева/справа от контактного разрыва
    double rho_sL = a_L / (S_L - S_C);
    double rho_sR = a_R / (S_R - S_C);

    // Удельная полная энергия слева/справа от контактного разрыва
    double E_sL = zL.E() + (S_C - u_L) * (S_C + P_L / a_L);
    double E_sR = zR.E() + (S_C - u_R) * (S_C + P_R / a_R);

    // Консервативный вектор слева/справа от контактного разрыва
    smf::QState Q_sL(rho_sL, {rho_sL * S_C, rho_sL * zL.v(), rho_sL * zL.w()}, rho_sL * E_sL);
    smf::QState Q_sR(rho_sR, {rho_sR * S_C, rho_sR * zR.v(), rho_sR * zR.w()}, rho_sR * E_sR);

    // Консервативный вектор слева/справа
    smf::QState Q_L(zL);
    smf::QState Q_R(zR);

    // Дифференциальный поток слева/справа
    smf::Flux F_L(zL);
    smf::Flux F_R(zR);

    smf::Flux F_sL = F_L.arr() - S_L * (Q_sL.arr() - Q_L.arr());
    smf::Flux F_sR = F_R.arr() - S_R * (Q_sR.arr() - Q_R.arr());

    return WaveConfig({.S_L = S_L, .S_C = S_C, .S_R = S_R,
                       .QsL = Q_sL, .FsL = F_sL,
                       .QsR = Q_sR, .FsR = F_sR});
}

HLLC::WaveConfig HLLC::wave_config(
        const phys::Eos& eosL, const smf::QState& Q_L, const smf::Flux& F_L,
        const phys::Eos& eosR, const smf::QState& Q_R, const smf::Flux& F_R) {

    Vector3d v_L = Q_L.momentum / Q_L.density;
    Vector3d v_R = Q_R.momentum / Q_R.density;
    double e_L = Q_L.energy / Q_L.density - 0.5 * v_L.squaredNorm();
    double e_R = Q_R.energy / Q_R.density - 0.5 * v_R.squaredNorm();

    // Скорость звука слева и справа
    double c_L = eosL.sound_speed_re(Q_L.density, e_L);
    double c_R = eosR.sound_speed_re(Q_R.density, e_R);

    // Оценки скоростей расходящихся волн
    double S_L = std::min({v_L.x() - c_L, v_R.x() - c_R, 0.0});
    double S_R = std::max({v_L.x() + c_L, v_R.x() + c_R, 0.0});

    // Поток через левую/правую волну
    smf::Flux a_L = F_L.arr() - S_L * Q_L.arr();
    smf::Flux a_R = F_R.arr() - S_R * Q_R.arr();

    // Скорость контактного разрыва
    double S_C = (a_R.momentum.x() - a_L.momentum.x()) / (a_R.mass - a_L.mass);

    // Плотность слева/справа от контактного разрыва
    double rho_sL = a_L.mass / (S_C - S_L);
    double rho_sR = a_R.mass / (S_C - S_R);

    // Удельная полная энергия слева/справа от контактного разрыва
    double E_sL = sqr(S_C) + (a_L.energy - S_C * a_L.momentum.x()) / a_L.mass;
    double E_sR = sqr(S_C) + (a_R.energy - S_C * a_R.momentum.x()) / a_R.mass;

    // Консервативный вектор слева/справа от контактного разрыва
    smf::QState Q_sL(rho_sL, {rho_sL * S_C, rho_sL * v_L.y(), rho_sL * v_L.z()}, rho_sL * E_sL);
    smf::QState Q_sR(rho_sR, {rho_sR * S_C, rho_sR * v_R.y(), rho_sR * v_R.z()}, rho_sR * E_sR);

    smf::Flux F_sL = a_L.arr() + S_L * Q_sL.arr();
    smf::Flux F_sR = a_R.arr() + S_R * Q_sR.arr();

#if 0
    double err1 = ((F_R .arr() - F_sR.arr()) - S_R * (Q_R .arr() - Q_sR.arr())).cwiseAbs().maxCoeff();
    double err2 = ((F_sR.arr() - F_sL.arr()) - S_C * (Q_sR.arr() - Q_sL.arr())).cwiseAbs().maxCoeff();
    double err3 = ((F_sL.arr() - F_L .arr()) - S_L * (Q_sL.arr() - Q_L .arr())).cwiseAbs().maxCoeff();
    double err4 = (F_L.arr() - F_R.arr() + S_R * (Q_R.arr() - Q_sR.arr()) +
            S_C * (Q_sR.arr() - Q_sL.arr()) + S_L * (Q_sL.arr() - Q_L.arr())).cwiseAbs().maxCoeff();

    if (std::max({err1, err2, err3, err4}) > 1.0e-6) {
        throw std::runtime_error("Rankine-Hugoniot conditions failed");
    }
#endif

    return WaveConfig({.S_L = S_L, .S_C = S_C, .S_R = S_R,
                       .QsL = Q_sL, .FsL = F_sL,
                       .QsR = Q_sR, .FsR = F_sR});
}

mmf::Flux HLLC::flux(const mmf::PState &zL, const mmf::PState &zR, const phys::MixturePT &mixture) const {
    return calc_flux(zL, zR, mixture);
}

mmf::Flux HLLC::calc_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::MixturePT &mixture) {
    using namespace mmf;

    const double &rho_L = zL.density;
    const double &rho_R = zR.density;
    const double &u_L = zL.velocity.x();
    const double &u_R = zR.velocity.x();
    const double &P_L = zL.pressure;
    const double &P_R = zR.pressure;

    // Скорость звука слева и справа
    Fractions a_L = zL.alpha();
    double c_L = mixture.sound_speed_rP(zL.density, zL.pressure, zL.mass_frac,
                                        {.T0 = zL.T(), .alpha = &a_L});
    Fractions a_R = zR.alpha();
    double c_R = mixture.sound_speed_rP(zR.density, zR.pressure, zR.mass_frac,
                                        {.T0 = zR.T(), .alpha = &a_R});

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
        PState z_sL(rho_sL, {S_C, zL.v(), zL.w()}, P_C, e_sL, NAN, zL.beta(), Fractions::NaN());

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
        PState z_sR(rho_sR, {S_C, zR.v(), zR.w()}, P_C, e_sR, NAN, zR.beta(), Fractions::NaN());

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

smf::Flux HLLC_LM::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux HLLC_LM::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    const double &rho_L = zL.density;
    const double &rho_R = zR.density;
    const double &u_L = zL.velocity.x();
    const double &u_R = zR.velocity.x();
    const double &P_L = zL.pressure;
    const double &P_R = zR.pressure;

    // Скорость звука слева и справа
    double c_L = eos.sound_speed_rP(rho_L, P_L);
    double c_R = eos.sound_speed_rP(rho_R, P_R);

    // Оценки скоростей расходящихся волн
    double S_L = std::min(u_L - c_L, u_R - c_R);
    double S_R = std::max(u_L + c_L, u_R + c_R);

    // Сверхзвуковое течение
    if (S_L >= 0.0) { return Flux(zL); }
    if (S_R <= 0.0) { return Flux(zR); }

    // Перенос массы через левую/правую волну
    double a_L = rho_L * (S_L - u_L);
    double a_R = rho_R * (S_R - u_R);

    // Скорость контактного разрыва
    double S_C = (P_L - P_R + a_R * u_R - a_L * u_L) / (a_R - a_L);

    // Плотность слева/справа от контактного разрыва
    double rho_sL = a_L / (S_L - S_C);
    double rho_sR = a_R / (S_R - S_C);

    // Удельная полная энергия слева/справа от контактного разрыва
    double E_sL = zL.E() + (S_C - u_L) * (S_C + P_L / a_L);
    double E_sR = zR.E() + (S_C - u_R) * (S_C + P_R / a_R);

    // Консервативный вектор слева/справа от контактного разрыва
    QState Q_sL(rho_sL, {rho_sL * S_C, rho_sL * zL.v(), rho_sL * zL.w()}, rho_sL * E_sL);
    QState Q_sR(rho_sR, {rho_sR * S_C, rho_sR * zR.v(), rho_sR * zR.w()}, rho_sR * E_sR);

    // Консервативный вектор слева/справа
    QState Q_L(zL);
    QState Q_R(zR);

    // Дифференциальный поток слева/справа
    Flux F_L(zL);
    Flux F_R(zR);

    // Индикатор, phi -> 0 при низком числе Маха (< 0.1)
    double Ma_limit = 0.1;
    double Ma_local = std::max(std::abs(u_L / c_L), std::abs(u_R / c_R));
    double xi = Ma_local / Ma_limit;
    double phi = xi > 1.0 ? 1.0 : std::sin(M_PI_2 * xi);
    //double phi = xi > 1.0 ? 1.0 : xi * (2.0 - xi);

    Flux F = 0.5 * (F_L.arr() + F_R.arr() +
            phi * S_L * (Q_sL.arr() - Q_L.arr()) +
            std::abs(S_C) * (Q_sL.arr() - Q_sR.arr()) +
            phi * S_R * (Q_sR.arr() - Q_R.arr())
    );

    if (F.arr().hasNaN()) {
        std::cerr << "HLLC_LM::calc_flux error\n";
        std::cerr << "  z_L: " << zL << "\n";
        std::cerr << "  z_R: " << zR << "\n";
        std::cerr << "  c_L: " << c_L << "; c_R: " << c_R << "\n";
        std::cerr << "  S_L: " << S_L << "; S_R: " << S_R << "\n";
        std::cerr << "  F_HLLC: " << F << "\n";
        throw std::runtime_error("HLLC_M::calc_flux error: bad value");
    }

    return F;
}

smf::Flux HLLC_M::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux HLLC_M::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    const double &rho_L = zL.density;
    const double &rho_R = zR.density;
    const double &u_L = zL.velocity.x();
    const double &u_R = zR.velocity.x();
    const double &P_L = zL.pressure;
    const double &P_R = zR.pressure;

    // Скорость звука слева и справа
    double c_L = eos.sound_speed_rP(rho_L, P_L);
    double c_R = eos.sound_speed_rP(rho_R, P_R);

    // Оценки скоростей расходящихся волн
    double S_L = std::min(u_L - c_L, u_R - c_R);
    double S_R = std::max(u_L + c_L, u_R + c_R);

    // Сверхзвуковое течение
    if (S_L >= 0.0) { return Flux(zL); }
    if (S_R <= 0.0) { return Flux(zR); }

    // Перенос массы через левую/правую волну
    double a_L = rho_L * (S_L - u_L);
    double a_R = rho_R * (S_R - u_R);

    // Скорость контактного разрыва
    double S_C = (P_L - P_R + a_R * u_R - a_L * u_L) / (a_R - a_L);

    // Касательные скорости также усредняются
    double V_C = (a_R * zR.v() - a_L * zL.v()) / (a_R - a_L);
    double W_C = (a_R * zR.w() - a_L * zL.w()) / (a_R - a_L);

    // Плотность слева/справа от контактного разрыва
    double rho_sL = a_L / (S_L - S_C);
    double rho_sR = a_R / (S_R - S_C);

    // Удельная полная энергия слева/справа от контактного разрыва
    double E_sL = zL.E() + (S_C - u_L) * (S_C + P_L / a_L);
    double E_sR = zR.E() + (S_C - u_R) * (S_C + P_R / a_R);

    // Квадраты тангенциальной составляющей скорости
    double v_L2 = sqr(zL.v()) + sqr(zL.w());
    double v_R2 = sqr(zR.v()) + sqr(zR.w());

    // Добавка для HLLC-M
    E_sL += 0.5 * ((a_R * v_R2 - a_L * v_L2) / (a_R - a_L) - v_L2);
    E_sR += 0.5 * ((a_R * v_R2 - a_L * v_L2) / (a_R - a_L) - v_R2);

    // Консервативный вектор слева/справа от контактного разрыва
    QState Q_sL(rho_sL, {rho_sL * S_C, rho_sL * V_C, rho_sL * W_C}, rho_sL * E_sL);
    QState Q_sR(rho_sR, {rho_sR * S_C, rho_sR * V_C, rho_sR * W_C}, rho_sR * E_sR);

    // Консервативный вектор слева/справа
    QState Q_L(zL);
    QState Q_R(zR);

    // Дифференциальный поток слева/справа
    Flux F_L(zL);
    Flux F_R(zR);

    // Симметричная формула HLLC
    Flux F = 0.5 * (F_L.arr() + F_R.arr() +
                    S_L * (Q_sL.arr() - Q_L.arr()) +
                    std::abs(S_C) * (Q_sL.arr() - Q_sR.arr()) +
                    S_R * (Q_sR.arr() - Q_R.arr())
    );

    if (F.arr().hasNaN()) {
        std::cerr << "HLLC_M::calc_flux error\n";
        std::cerr << "  z_L: " << zL << "\n";
        std::cerr << "  z_R: " << zR << "\n";
        std::cerr << "  c_L: " << c_L << "; c_R: " << c_R << "\n";
        std::cerr << "  S_L: " << S_L << "; S_R: " << S_R << "\n";
        std::cerr << "  F_HLLC: " << F << "\n";
        throw std::runtime_error("HLLC_M::calc_flux error: bad value");
    }

    return F;
}

}