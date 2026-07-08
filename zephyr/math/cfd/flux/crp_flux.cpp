#include <iostream>
#include <zephyr/math/cfd/flux/crp_flux.h>
#include <zephyr/math/cfd/flux/hll.h>
#include <zephyr/math/cfd/flux/hllc.h>


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
        if (ts < t) {
            return std::numeric_limits<double>::infinity();
        } else {
            return ts;
        }
    }

    // Пересечение характеристик, находятся только пересечения во времени дальше,
    // чем время обеих точек исходных характеристик. В обратном случае считается,
    // что характеристики пересекаются на бесконечности t = +inf.
    Point cross(const Char &c) const {
        constexpr double inf = std::numeric_limits<double>::infinity();

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

//#define USE_HLL_FLUX
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

    if (S_0L > 0.0) {
        throw std::runtime_error("Supersonic HLL flux #1");
    }

    // Первый поток через грань
    const mmf::Flux& F1 = F_s1;

    // Характеристика начального контакта
    Char C_0 = {.x = -delta, .t = 0.0, .S = zLA.velocity.x()};

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
    Char C_1C = {.x = O1.x, .t = O1.t, .S = S_1C};
    Char C_1R = {.x = O1.x, .t = O1.t, .S = S_1R};

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

// Очень простая функция, выбирает решение в нужной области в зависимости от S
inline std::tuple<mmf::QState, mmf::Flux> choose(
    double S, double S_L, double S_C, double S_R,
    const mmf::QState& Q_L, const mmf::QState& Q_sL, const mmf::QState& Q_sR, const mmf::QState& Q_R,
    const mmf::Flux&   F_L, const mmf::Flux&   F_sL, const mmf::Flux&   F_sR, const mmf::Flux&   F_R) {
    if (S < S_C) {
        if (S < S_L)
            return {Q_L,  F_L};
        else
            return {Q_sL, F_sL};
    }
    else {
        if (S < S_R)
            return {Q_sR, F_sR};
        else
            return {Q_R,  F_R};
    }
}

// Классическая задача для CRP
mmf::Flux CrpFlux::classic(const mmf::PState& zLA, const mmf::PState& zLB, const mmf::PState& zRB,
                           const MixturePT& mixture, double delta, double dt) {

#if 0 // MRV VERSION

    mmf::QState qLB(zLB);
    mmf::Flux   fLB(zLB);
    mmf::QState Q_R(zRB);
    mmf::Flux   F_R(zRB);

    // Характеристика начального контакта
    Char C_0 = {.x = -delta, .t = 0.0, .S = zLA.vx()};

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
        Char C_1C = {.x = O1.x, .t = O1.t, .S = S_1C};
        Char C_1R = {.x = O1.x, .t = O1.t, .S = S_1R};

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
        Char C_1R = {.x = O1.x, .t = O1.t, .S = S_1R};

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

#else // ZPP VERSION

    // Начало координат
    double x0{0.0}, t0{0.0};

    // Первый распад
    mmf::QState qLB(zLB);
    mmf::Flux   fLB(zLB);
    mmf::QState Q_R(zRB);
    mmf::Flux   F_R(zRB);
    auto[S_0L, S_0C, S_0R, Q_s0L, F_s0L, Q_s0R, F_s0R] = HLLC::wave_config(mixture, qLB, fLB, Q_R, F_R);

    // точка O1 = (x1, t1) - взаимодействие C_0 и C_0L
    Char C_0  = {.x = -delta, .t = 0.0, .S = zLA.vx()};
    Char C_0L = {.x = 0.0, .t = 0.0, .S = S_0L};
    auto [x1, t1]  = C_0.cross(C_0L);

    // Сверхзвук, точка O1 в правой ячейке
    if( x1 >= 0.0 ) {
        double tau_0 = C_0.edge_t();
        if( tau_0 >= dt ) {
            return fLB;
        }

        mmf::Flux G_O0_O1 = fLB.arr()*( t1 - t0 ) - qLB.arr()*( x1 - x0 );

        // Второй распад
        mmf::QState qLA(zLA);
        mmf::Flux   fLA(zLA);

        auto F = fLA;
        auto Q = qLA;

        if( dt > t1 ) {
            auto [S_1L, S_1C, S_1R, Q_s1L, F_s1L, Q_s1R, F_s1R] = HLLC::wave_config(mixture, qLA, fLA, Q_s0L, F_s0L);

            auto lb = -x1 / (dt - t1);
            std::tie(Q, F) = choose(lb, S_1L, S_1C, S_1R,
                                    qLA, Q_s1L, Q_s1R, Q_s0L,
                                    fLA, F_s1L, F_s1R, F_s0L);
        }

        mmf::Flux G_O1_Odt = F.arr()*( dt - t1 ) - Q.arr()*( x0 - x1 );

        return (G_O0_O1.arr() + G_O1_Odt.arr()) / dt;
    }

    // mix2clear
    if( S_0C >= 0.0 ) {
        // Нет взаимодействия S_0L за dt
        if( t1 >= dt ) {
            return F_s0L;
        }

        // Второй распад
        mmf::QState Q_L(zLA);
        mmf::Flux   F_L(zLA);
        auto [S_1L, S_1C, S_1R, Q_s1L, F_s1L, Q_s1R, F_s1R] = HLLC::wave_config(mixture, Q_L, F_L, Q_s0L, F_s0L);

        // Формально захватывается случаем Odt O1 O
        Char C_1R = {.x = x1, .t = t1, .S = S_1R};
        double tau_1 = C_1R.edge_t();
        if( tau_1 >= dt ) {
            return F_s0L;
        }

        // Версия без разделения потоков материалов A/B. Треугольник: Odt O1 O0.
        // Я проверил версию с разделением потоков, математически формулы эквивалентны,
        // но по какой-то причине она оказывается менее устойчивой?
        mmf::Flux G_O0_O1 = F_s0L.arr()*( t1 - t0 ) - Q_s0L.arr()*( x1 - x0 );

        double lb = x1 / ( t1 - dt );
        auto[Q, F] = choose(lb, S_1L, S_1C, S_1R,
                            Q_L, Q_s1L, Q_s1R, Q_s0L,
                            F_L, F_s1L, F_s1R, F_s0L);

        mmf::Flux G_01_Odt = F.arr()*( dt - t1 ) - Q.arr()*( x0 - x1 );

        return ( G_O0_O1.arr() + G_01_Odt.arr() ) / dt;
    }
    // clear2mix
    else {
        // Нет взаимодействия S_0L за dt
        if( t1 >= dt ) {
            return F_s0R;
        }

        // Второй распад в точке O1
        mmf::QState Q_L(zLA);
        mmf::Flux   F_L(zLA);
        auto[S_1L, S_1C, S_1R, Q_s1L, F_s1L, Q_s1R, F_s1R] = HLLC::wave_config(mixture, Q_L, F_L, Q_s0L, F_s0L);

        // Точка O2 = (x2, t2) - взаимодействие C_0C и C_1R
        Char C_0C = {.x = 0.0,  .t = 0.0,  .S = S_0C};
        Char C_1R = {.x = x1, .t = t1, .S = S_1R};
        auto [x2, t2]  = C_1R.cross(C_0C);

        // нет взаимодействия C_1R C_0C за dt
        if( t2 > dt ) {
            return F_s0R;
        }

        // Редкий случай при скорости S_0С около нуля
        if( x2 > 0.0 ){
            // При скорости S_0C около нуля поток непрерывен, то есть F_s0R ~ F_s0L,
            // поэтому вообще говоря не принципиально какой использовать.
            // Но почему распады далее не рассматриваются? Свести к предыдущим if-ам
            return F_s0R;
        }

        // Третий распад в точке O2
        auto[S_2L, S_2C, S_2R, Q_s2L, F_s2L, Q_s2R, F_s2R] = HLLC::wave_config(mixture, Q_s1R, F_s1R, Q_s0R, F_s0R);

        // Точка O3 = (x3, t3) - взаимодействие C_2L и C_1C
        Char C_2L = {.x = x2, .t = t2, .S = S_2L};
        Char C_1C = {.x = x1, .t = t1, .S = S_1C};
        auto [x3, t3] = C_1C.cross(C_2L);

        //case 2 search
        if( x3 >= 0.0 ) {
            throw std::runtime_error("Flux::classic: very strange case");
        }

        // Четвертый распад в точке 03
        auto[S_3L, S_3C, S_3R, Q_s3L, F_s3L, Q_s3R, F_s3R] = HLLC::wave_config(mixture, Q_s1L, F_s1L, Q_s2L, F_s2L);

        Char C_3R = {.x = x3, .t = t3, .S = S_3R};
        double tau_3 = C_3R.edge_t();

        // Вещество не перетекает
        if( dt <= tau_3 ) {
            // Формально этот случай захватывается формулами ниже
            // внесен в этот иф, тк может быть t1 > t4 ????
            Char C_2R = {.x = x2, .t = t2, .S = S_2R};
            double tau_2 = C_2R.edge_t();
            if( dt <= tau_2 ) {
                return F_s0R;
            }

            // Треугольник: O0 O2 Odt
            mmf::Flux G_O0_O2 = F_s0R.arr()*( t2 - t0 ) - Q_s0R.arr()*( x2 - x0 );

            // Гарантируется x2 < 0, t2 < dt, т.е. lb > 0
            double lb = -x2 / (dt - t2);
            auto[Q, F] = choose(lb, S_2L, S_2C, S_2R,
                                Q_s1R, Q_s2L, Q_s2R, Q_s0R,
                                F_s1R, F_s2L, F_s2R, F_s0R);

            mmf::Flux G_O2_Odt = F.arr()*( dt - t2 ) - Q.arr()*( x0 - x2 );

            return ( G_O0_O2.arr() + G_O2_Odt.arr() ) / dt;
        }
        else {
            // Полигон: O0 O2 O3 Odt
            mmf::Flux G_O0_O2 = F_s0R.arr()*( t2 - t0 ) - Q_s0R.arr()*( x2 - x0 );
            mmf::Flux G_O2_O3 = F_s2L.arr()*( t3 - t2 ) - Q_s2L.arr()*( x3 - x2 );

            double lb = -x3/(dt - t3);
            auto[Q, F] = choose(lb, S_3L, S_3C, S_3R,
                                Q_s1L, Q_s3L, Q_s3R, Q_s2L,
                                F_s1L, F_s3L, F_s3R, F_s2L);

            mmf::Flux G_O3_Odt = F.arr()*( dt - t3 ) - Q.arr()*( x0 - x3 );

            return ( G_O0_O2.arr() + G_O2_O3.arr() + G_O3_Odt.arr() ) / dt;
        }
    }
#endif
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