#include <iostream>
#include <boost/function/function_base.hpp>
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

#define USE_HLL_FLUX
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

    //ПЕРВЫЙ распад
     mmf::QState qLB(zLB);
     mmf::Flux   fLB(zLB);
     mmf::QState Q_R(zRB);
     mmf::Flux   F_R(zRB);
     auto[S_0L, S_0C, S_0R, Q_s0L, F_s0L, Q_s0R, F_s0R] = HLLC::wave_config(mixture, qLB, fLB, Q_R, F_R);


     //точка O1, OO, Odt
     Char C_0  = {.x = -delta, .t = 0.0, .S = zLA.vx()};
     Char C_0L = {.x = 0.0, .t = 0.0, .S = S_0L};
     Point O1  = C_0.cross(C_0L);
     Point O0( 0.0, 0.0);
     Point Odt(0.0, dt );

     //mix2clear
     if( S_0C >= 0.0 ) {

         // Нет взаимодействия S_0L за dt
         if( dt <= O1.t ) {
             return F_s0L;
         }

         // ВТОРОЙ распад
         mmf::QState Q_L(zLA);
         mmf::Flux   F_L(zLA);
         auto [S_1L, S_1C, S_1R, Q_s1L, F_s1L, Q_s1R, F_s1R] = HLLC::wave_config(mixture, Q_L, F_L, Q_s0L, F_s0L);

         //time t1
         Char C_1R = {.x = O1.x, .t = O1.t, .S = S_1R};
         double t1 = C_1R.edge_t();

         //формально захватывается случаем Odt O1 O
         if( dt <= t1 ) {
             return F_s0L;
         }

         //version without distinguishing A/B material fluxes
         //triangle: Odt O1 O
         //full flux Ftot
         double lb = -O1.x/( dt - O1.t );
         if( lb > S_1R ) {
             throw std::runtime_error("Flux::classic: larger than S_1L");
         }
         auto F =  lb < S_1L ? F_L : ( lb < S_1C ? F_s1L : F_s1R );
         auto Q =  lb < S_1L ? Q_L : ( lb < S_1C ? Q_s1L : Q_s1R );
         //auto F =                    lb < S_1C ? F_s1L : ( lb < S_1L ? F_s1R : F_s0L ); //избыточный в силу if else
         //auto Q =                    lb < S_1C ? Q_s1L : ( lb < S_1L ? Q_s1R : Q_s0L );
         mmf::Flux G_OdtO1 = F.arr()*( O1.t  - Odt.t ) - Q.arr()*( O1.x - Odt.x );

         F = F_s0L;
         Q = Q_s0L;
         mmf::Flux G_O1O = F.arr()*( O0.t - O1.t ) - Q.arr()*( O0.x - O1.x );

         mmf::Flux Ftot1 = ( G_OdtO1.arr() + G_O1O.arr() )/(0.0 - dt);
         //return Ftot1;


         //version to distinguish A/B material fluxes
         //triangles:  Ott O1 O  and  Odt O1 Ott

         //flux of B of Ftot
         Char C_1C = {.x = O1.x, .t = O1.t, .S = S_1C};
         double t2 = C_1C.edge_t();
         double tt = std::min(dt, t2);
         double w1 = tt/dt;

         Point Ott(0.0, tt);
         lb = -O1.x/( tt - O1.t );
         F = F_s1R;
         Q = Q_s1R;
         mmf::Flux G_OttO1 = F.arr()*( O1.t  - Ott.t ) - Q.arr()*( O1.x - Ott.x );

         F = F_s0L;
         Q = Q_s0L;
         G_O1O = F.arr()*( O0.t - O1.t ) - Q.arr()*( O0.x - O1.x );

         mmf::Flux F_B = ( G_OttO1.arr() + G_O1O.arr() )/( 0.0 - tt );

         if( dt < t2) {
             return F_B;
         }


         //flux of A of Ftot
         Point Ot2(0.0, t2);
         lb = -O1.x/( dt - O1.t );
         F = lb < S_1L  ?  F_L  :  F_s1L;
         Q = lb < S_1L  ?  Q_L  :  Q_s1L;
         G_OdtO1 = F.arr()*( O1.t  - Ott.t ) - Q.arr()*( O1.x - Ott.x );

         F = F_s1L;
         Q = Q_s1L;
         mmf::Flux G_O1Ot2 = F.arr()*( Ot2.t  - O1.t ) - Q.arr()*( Ot2.x - O1.x );

         mmf::Flux F_A = ( G_OdtO1.arr() + G_O1Ot2.arr() )/( t2 - dt );

         mmf::Flux Ftot = w1*F_B.arr() + (1.0 - w1)*F_A.arr();
         return Ftot;

     }//clear2mix
     else{

         // Нет взаимодействия S_0L за dt
         if( O1.t >= dt ) {
             return F_s0R;
         }


         //ВТОРОЙ распад в точке O1
         mmf::QState Q_L(zLA);
         mmf::Flux   F_L(zLA);
         auto[S_1L, S_1C, S_1R, Q_s1L, F_s1L, Q_s1R, F_s1R] = HLLC::wave_config(mixture, Q_L, F_L, Q_s0L, F_s0L);

         //точка O2 - взаимодействие C_1R C_0C
         Char C_0C = {.x = 0.0,  .t = 0.0,  .S = S_0C};
         Char C_1R = {.x = O1.x, .t = O1.t, .S = S_1R};
         Point O2  = C_1R.cross(C_0C);

         // нет взаимодействия C_1R C_0C за dt
         // или редкая хуета при малых C_0C.S; O2.x > 0.0
         if( O2.t > dt || O2.x > 0.0 ) {
             return F_s0R;
         }


         //ТРЕТИЙ распад в точке O2
         auto[S_2L, S_2C, S_2R, Q_s2L, F_s2L, Q_s2R, F_s2R] = HLLC::wave_config(mixture, Q_s1R, F_s1R, Q_s0R, F_s0R);

         //точка O3
         Char C_2L = {.x = O2.x, .t = O2.t, .S = S_2L};
         Char C_1C = {.x = O1.x, .t = O1.t, .S = S_1C};
         Point O3 = C_1C.cross(C_2L);

         //case 2 search
         if( O3.x >= 0.0 ) {
             throw std::runtime_error("Flux::classic: c2m;  very strange case");
         }

         //ЧЕТВЕРТНЫЙ распад в точке 03
         auto[S_3L, S_3C, S_3R, Q_s3L, F_s3L, Q_s3R, F_s3R] = HLLC::wave_config(mixture, Q_s1L, F_s1L, Q_s2L, F_s2L);

         //some times
         Char C_2C = {.x = O2.x, .t = O2.t, .S = S_2C};
         double t3 = C_2C.edge_t();
         Char C_3C = {.x = O3.x, .t = O3.t, .S = S_3C};
         double t5 = C_3C.edge_t();

         Char C_2R = {.x = O2.x, .t = O2.t, .S = S_2R};
         Char C_3R = {.x = O3.x, .t = O3.t, .S = S_3R};
         double t1 = C_2R.edge_t();
         double t4 = C_3R.edge_t();


         if( dt <= t4 ) {
         //by triangle Odt O2 O

             //формально этот случай захватывается формулами ниже
             //внесен в этот иф, тк может быть t1 > t4
             if( dt <= t1 ) {
                 return F_s0R;
             }

             double lb = -O2.x/( dt - O2.t );
             //auto F = lb < S_2C ? F_s2L :               F_s2R;
             //auto Q = lb < S_2C ? Q_s2L :               Q_s2R;
             auto   F = lb < S_2C ? F_s2L : ( lb < S_2R ? F_s2R : F_s0R );
             auto   Q = lb < S_2C ? Q_s2L : ( lb < S_2R ? Q_s2R : Q_s0R );
             mmf::Flux G_OdtO2 = F.arr()*( O2.t  - Odt.t ) - Q.arr()*( O2.x - Odt.x);

             F = F_s0R;
             Q = Q_s0R;
             mmf::Flux G_O2O0   = F.arr()*( O0.t - O2.t )  - Q.arr()*( O0.x - O2.x );

             mmf::Flux Ftot = ( G_OdtO2.arr() + G_O2O0.arr() )/(0.0 - dt);
             return Ftot;

         }else{
         //by poly Odt O3 O2 O

             double lb = -O3.x/(dt - O3.t);

             //auto F =                     lb < S_3C ? F_s3L : F_s3R;
             //auto Q =                     lb < S_3C ? Q_s3L : Q_s3R;
             auto F = lb < S_3L ? F_s1L : ( lb < S_3C ? F_s3L : F_s3R );
             auto Q = lb < S_3L ? Q_s1L : ( lb < S_3C ? Q_s3L : Q_s3R );
             mmf::Flux G_OdtO3 = F.arr()*( O3.t  - Odt.t ) - Q.arr()*( O3.x - Odt.x);

             F = F_s2L;
             Q = Q_s2L;
             mmf::Flux G_O3O2  = F.arr()*( O2.t  - O3.t )  - Q.arr()*( O2.x - O3.x );

             F = F_s0R;
             Q = Q_s0R;
             mmf::Flux G_O2O0  = F.arr()*( O0.t - O2.t )   - Q.arr()*( O0.x - O2.x );

             mmf::Flux Ftot = ( G_OdtO3.arr() + G_O3O2.arr() + G_O2O0.arr() )/(0.0 - dt);
             return Ftot;
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