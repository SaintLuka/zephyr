#include <zephyr/math/cfd/flux/crp_flux.h>
#include <zephyr/math/cfd/flux/hll.h>
#include <zephyr/math/cfd/flux/hllc.h>

//#define USE_HLL_FLUX

namespace zephyr::math {

using zephyr::phys::MixturePT;

namespace {

// Точка в плоскости (x, t)
struct Point {
    double x, t;
};

// Характеристика (x, t) -- точка, S -- скорость
struct Char {
    double x, t, S;

    // Пересечение характеристики с прямой x = 0, время t.
    double edge_t() const {
        double ts = (S * t - x) / S;
        if (ts <= t) {
            return std::numeric_limits<double>::infinity();
        } else {
            return ts;
        }
    }

    // Пересечение характеристик, находятся только пересечения во времени дальше,
    // чем время обеих точек исходных характеристик. В обратном случае считается,
    // что характеристики пересекаются на бесконечности t = +inf.
    Point cross(const Char &c) const {
        const double inf = std::numeric_limits<double>::infinity();

        if (S == c.S) {
            return {.x = 0.5 * (x + c.x), .t = inf};
        }

        Point p = {.x = NAN, .t = NAN};
        p.t = (c.x - x + S * t - c.S * c.t) / (S - c.S);

        if (p.t <= std::max(t, c.t)) {
            return {.x = 0.5 * (x + c.x), .t = inf};
        }

        p.x = 0.5 * (x + c.x + S * (p.t - t) + c.S * (p.t - c.t));

        return p;
    }
};

}

#ifdef USE_HLL_FLUX

// Классическая задача для CRP
mmf::Flux CrpFlux::classic(const mmf::PState& zLA, const mmf::PState& zLB, const mmf::PState& zRB,
                           const MixturePT& mixture, double delta, double dt) {
    mmf::QState qLB(zLB);
    mmf::Flux   fLB(zLB);
    mmf::QState Q_R(zRB);
    mmf::Flux   F_R(zRB);

    // Характеристики из HLL
    auto[S_0L, S_0R, Q_s1, F_s1] = HLL::wave_config(mixture, qLB, fLB, Q_R, F_R);
    Char C_0L = {.x = 0.0, .t = 0.0, .S = S_0L};
    Char C_0R = {.x = 0.0, .t = 0.0, .S = S_0R};

    // Первый поток через грань
    const mmf::Flux& F1 = F_s1;

    // Характеристика начального контакта
    Char C_0 = {.x = -delta, .t = 0.0, .S = zLA.u()};

    // Первое взаимодействие
    Point O1 = C_0.cross(C_0L);

    // Нет взаимодействия за dt
    if (O1.t >= dt) {
        return F1;
    }

    mmf::QState Q_L(zLA);
    mmf::Flux   F_L(zLA);

    // Пара характеристик из HLLC
    auto[S_1L, S_1C, S_1R, Q_s2L, F_s2L, Q_s2R, F_s2R] = HLLC::wave_config(mixture, Q_L, F_L, Q_s1, F_s1);
    Char C_1C = {.x = O1.x, .t = O1.t, S_1C};
    Char C_1R = {.x = O1.x, .t = O1.t, S_1R};

    double tau1 = C_1R.edge_t() / dt;
    if (tau1 >= 1.0) {
        return F1;
    }

    // Второй поток
    const mmf::Flux& F2 = F_s2R;

    // Второе взаимодействие
    //Point O2 = C_1R & C_0R;

    double tau2 = C_1C.edge_t() / dt;


    if (tau2 >= 1.0) {
        return tau1 * F1.arr() + (1.0 - tau1) * F2.arr();
    }

    const mmf::Flux& F3 = F_s2L;

    return tau1 * F1.arr() + (tau2 - tau1) * F2.arr() + (1.0 - tau2) * F3.arr();
}

#else

// Классическая задача для CRP
mmf::Flux CrpFlux::classic(const mmf::PState& zLA, const mmf::PState& zLB, const mmf::PState& zRB,
                           const MixturePT& mixture, double delta, double dt) {
    mmf::QState qLB(zLB);
    mmf::Flux   fLB(zLB);
    mmf::QState Q_R(zRB);
    mmf::Flux   F_R(zRB);

    // Характеристика начального контакта
    Char C_0 = {.x = -delta, .t = 0.0, .S = zLA.u()};

    // Характеристики из точки O
    auto[S_0L, S_0C, S_0R, Q_s0L, F_s0L, Q_s0R, F_s0R] = HLLC::wave_config(mixture, qLB, fLB, Q_R, F_R);
    Char C_0L = {.x = 0.0, .t = 0.0, .S = S_0L};

    // Первое взаимодействие
    Point O1 = C_0.cross(C_0L);

    if (S_0C >= 0.0) {
        // Положительная скорость в веществе 'B'.

        // Первый поток через грань
        const mmf::Flux& F0 = F_s0L;

        // Нет взаимодействия за dt
        if (O1.t >= dt) {
            return F0;
        }

        mmf::QState Q_L(zLA);
        mmf::Flux   F_L(zLA);

        // Характеристики из точки O1
        auto[S_1L, S_1C, S_1R, Q_s1L, F_s1L, Q_s1R, F_s1R] = HLLC::wave_config(mixture, Q_L, F_L, Q_s0L, F_s0L);
        Char C_1C = {.x = O1.x, .t = O1.t, S_1C};
        Char C_1R = {.x = O1.x, .t = O1.t, S_1R};

        double tau1 = C_1R.edge_t() / dt;
        if (tau1 >= 1.0) {
            return F0;
        }

        // Второй поток
        const mmf::Flux& F1 = F_s1R;

        double tau2 = C_1C.edge_t() / dt;

        if (tau2 >= 1.0) {
            return tau1 * F0.arr() + (1.0 - tau1) * F1.arr();
        }

        const mmf::Flux& F2 = F_s1L;

        return tau1 * F0.arr() + (tau2 - tau1) * F1.arr() + (1.0 - tau2) * F2.arr();
    }
    else {
        // Отрицательная скорость в веществе 'B'.

        // Первый поток через грань
        const mmf::Flux& F0 = F_s0R;

        // Нет взаимодействия за dt
        if (O1.t >= dt) {
            return F0;
        }

        mmf::QState Q_L(zLA);
        mmf::Flux   F_L(zLA);

        // Контактный разрыв в веществе B.
        Char C_0C = {.x = 0.0, .t = 0.0, .S = S_0C};

        // Характеристики из точки O1
        auto[S_1L, S_1C, S_1R, Q_s1L, F_s1L, Q_s1R, F_s1R] = HLLC::wave_config(mixture, Q_L, F_L, Q_s0L, F_s0L);
        Char C_1R = {.x = O1.x, .t = O1.t, S_1R};

        Point O2 = C_1R.cross(C_0C);

        if (O2.t > dt) {
            return F0;
        }

        // Характеристики из точки O2
        auto[S_2L, S_2C, S_2R, Q_s2L, F_s2L, Q_s2R, F_s2R] = HLLC::wave_config(mixture, Q_s1R, F_s1R, Q_s0R, F_s0R);
        Char C_2C = {.x = O2.x, .t = O2.t, .S = S_2C};
        Char C_2R = {.x = O2.x, .t = O2.t, .S = S_2R};

        double tau1 = C_2R.edge_t() / dt;
        if (tau1 >= 1.0) {
            return F0;
        }

        const mmf::Flux& F1 = F_s2R;

        double tau2 = C_2C.edge_t() / dt;
        if (tau2 >= 1.0) {
            return tau1 * F0.arr() + (1.0 - tau1) * F1.arr();
        }

        // Экзотический случай
        const mmf::Flux& F2 = F_s2L;

        mmf::Flux res = tau1 * F0.arr() + (tau2 - tau1) * F1.arr() + (1.0 - tau2) * F2.arr();

        // Нужно проверить, что всё вещество не вытекает из левой
        return res;
    }
}

#endif

mmf::Flux CrpFlux::inverse(const mmf::PState& zLA, const mmf::PState& zRA, const mmf::PState& zRB,
                           const MixturePT& mixture, double delta, double dt) {
    auto& zLA_i = const_cast<mmf::PState&>(zLA);
    auto& zRA_i = const_cast<mmf::PState&>(zRA);
    auto& zRB_i = const_cast<mmf::PState&>(zRB);

    zLA_i.inverse();
    zRA_i.inverse();
    zRB_i.inverse();

    auto flux = CrpFlux::classic(zRB_i, zRA_i, zLA_i, mixture, delta, dt);
    flux.inverse();

    zLA_i.inverse();
    zRA_i.inverse();
    zRB_i.inverse();

    return flux;
}

} // namespace zephyr::math