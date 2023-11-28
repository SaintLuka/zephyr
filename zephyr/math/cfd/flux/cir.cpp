#include <iostream>

#include <zephyr/geom/vector.h>
#include <zephyr/math/cfd/flux/cir.h>

namespace zephyr::math {

using geom::Vector5d;
using geom::Vector6d;
using geom::Matrix5d;

inline double sqr(double x) {
    return x * x;
}

smf::Flux CIR1::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux CIR1::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    // Состояние на грани
    PState zE = 0.5 * (zL.vec() + zR.vec());

    // Согласуем состояние на грани
    zE.energy = eos.energy_rp(zE.density, zE.pressure);

    auto P = eos.pressure_re(zE.density, zE.energy, {.deriv = true});

    double rho = zE.density;
    double u = zE.velocity.x();
    double v = zE.velocity.y();
    double w = zE.velocity.z();
    double p = zE.pressure;
    double eps = zE.energy;

    double u2 = sqr(u);
    double v2 = sqr(v);
    double w2 = sqr(w);

    double q2 = u2 + v2 + w2; // квадрат модуля скорости

    double E = rho * (eps + 0.5 * q2);

    double c = eos.sound_speed_rp(rho, p);
    double c2 = sqr(c);

    double h = (E + p) / rho; // энтальпия
    double b = P.dE / rho;
    double T = b * (q2 - h) + c2;

    // @formatter:off

    // Вектор с собственными значениями
    Vector5d vec_L;
    vec_L << u - c, u, u, u, u + c;

    // Матрица с модулями собственных значений на диагонали
    DiagMatrix5d abs_L = vec_L.cwiseAbs().asDiagonal();

    // Матрица слева в спектральном разложении якобиана.
    // Правые собственные вектора якобиана по столбцам
    Matrix5d OmR;
    OmR << 1          ,  0  ,  0  ,  1           ,  1          ,
           u - c      ,  0  ,  0  ,  u           ,  u + c      ,
           v          ,  1  ,  0  ,  v           ,  v          ,
           w          ,  0  ,  1  ,  w           ,  w          ,
           h - u * c  ,  v  ,  w  ,  h - c2 / b  ,  h + u * c  ;

    // Матрица справа в спектральном разложении якобиана.
    // Левые собственные вектора якобиана по строкам
    Matrix5d OmL;
    OmL << T + u * c         , -b * u - c  , -b * v      , -b * w      ,  b      ,
          -2 * v * c2        ,  0          ,  2 * c2     ,  0          ,  0      ,
          -2 * w * c2        ,  0          ,  0          ,  2 * c2     ,  0      ,
           2 * b * (h - q2)  ,  2 * b * u  ,  2 * b * v  ,  2 * b * w  , -2 * b  ,
           T - u * c         , -b * u + c  , -b * v      , -b * w      ,  b      ;
    OmL /= (2.0 * c2);

#if 0
    // Диагональ с собственными значениями
    DiagMatrix5d L = vec_L.asDiagonal();

    // Якобиан на грани
    Matrix5d A;
    A <<  0            ,  1              ,  0          ,  0          ,  0           ,
          T - u2       ,  (2 - b) * u    , -b * v      , -b * w      ,  b           ,
         -u * v        ,  v              ,  u          ,  0          ,  0           ,
         -u * w        ,  w              ,  0          ,  u          ,  0           ,
          u * (T - h)  ,  h - b * u2     , -b * u * v  , -b * u * w  ,  (1 + b) * u ;

    double norm1 = (OmR * L * OmL - A).norm();
    double norm2 = (OmL*OmR - Matrix5d::Identity()).norm();
    double norm3 = (OmR*OmL - Matrix5d::Identity()).norm();

    double err = 1.0e-10 * A.norm();

    // Проверяем, что OmL и OmR взаимно-обратные
    if (norm1 > err || norm2 > err || norm3 > err) {
        std::cerr << "Ошибка в декомпозиции\n";

        std::cerr << "||OmR * L * OmL - A||: " << norm1 << "\n";
        std::cerr << "||OmL * OmR - I||: " << norm2 << "\n";
        std::cerr << "||OmR * OmL - I||: " << norm3 << "\n";
        throw std::runtime_error("CIR: error");
    }
#endif

    // @formatter:on

    Matrix5d abs_A = OmR * abs_L * OmL;

    QState qL(zL); // Консервативный вектор слева
    QState qR(zR); // Консервативный вектор справа

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа

    Flux res = 0.5 * (fL.vec() + fR.vec()) + 0.5 * abs_A * (qL.vec() - qR.vec());

    return res;
}

smf::Flux CIR2::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux CIR2::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    // Состояние на грани
    PState zE = 0.5 * (zL.vec() + zR.vec());

    // Согласуем состояние на грани
    zE.energy = eos.energy_rp(zE.density, zE.pressure);

    double rho = zE.density;
    double u = zE.velocity.x();
    double v = zE.velocity.y();
    double w = zE.velocity.z();
    double p = zE.pressure;
    double eps = zE.energy;

    double u2 = sqr(u);
    double v2 = sqr(v);
    double w2 = sqr(w);

    double q2 = u2 + v2 + w2;

    double E = rho * (eps + 0.5 * q2);

    double c = eos.sound_speed_rp(rho, p);
    double c2 = sqr(c);

    double h = (E + p) / rho;

    // @formatter:off

    // Вектор с собственными значениями
    Vector6d vec_L;
    vec_L << u - c, u, u, u, u + c, u;

    // Матрица с модулями собственных значений на диагонали
    DiagMatrix6d abs_L = vec_L.cwiseAbs().asDiagonal();

    // Матрица слева в спектральном разложении якобиана.
    // Правые собственные вектора якобиана по столбцам
    Matrix<double, 5, 6> OmR;
    OmR << 1          ,  0  ,  0  ,  0  ,  1          ,  1  ,
           u - c      ,  0  ,  0  ,  0  ,  u + c      ,  u  ,
           v          ,  1  ,  0  ,  0  ,  v          ,  0  ,
           w          ,  0  ,  1  ,  0  ,  w          ,  0  ,
           h - u * c  ,  0  ,  0  ,  1  ,  h + u * c  ,  0  ;

    // Матрица справа в спектральном разложении якобиана.
    // Левые собственные вектора якобиана по строкам
    Matrix6d OmL;
    OmL <<  u / (2*c)  , -1 / (2*c)  ,  0  ,  0  ,  0  ,  1 / (2*c2)  ,
            0          ,  0          ,  1  ,  0  ,  0  , -v / c2      ,
            0          ,  0          ,  0  ,  1  ,  0  , -w / c2      ,
            u2         , -u          ,  0  ,  0  ,  1  , -h / c2      ,
           -u / (2*c)  ,  1 / (2*c)  ,  0  ,  0  ,  0  ,  1 / (2*c2)  ,
            1          ,  0          ,  0  ,  0  ,  0  , -1 / c2      ;

#if 1
    // Диагональ с собственными значениями
    DiagMatrix6d L = vec_L.asDiagonal();

    // Якобиан на грани
    Matrix<double, 5, 6> A;
    A <<  0      ,  1      ,  0  ,  0  ,  0  ,  0  ,
         -u2     ,  2 * u  ,  0  ,  0  ,  0  ,  1  ,
         -u * v  ,  v      ,  u  ,  0  ,  0  ,  0  ,
         -u * w  ,  w      ,  0  ,  u  ,  0  ,  0  ,
         -u * h  ,  h      ,  0  ,  0  ,  u  ,  u  ;

    double norm1 = (OmR * L * OmL - A).norm();
    Matrix<double, 5, 6> I2 = Matrix<double, 5, 6>::Identity();
    double norm2 = (OmR*OmL - I2).norm();

    double err = 1.0e-10 * A.norm();

    // Проверяем, что OmL и OmR взаимно-обратные
    if (norm1 > err || norm2 > err) {
        std::cerr << "Ошибка в декомпозиции\n";

        std::cerr << "||OmR * L * OmL - A||: " << norm1 << "\n";
        std::cerr << "||OmR * OmL - I||: " << norm2 << "\n";
        throw std::runtime_error("CIR: error");
    }
#endif

    // @formatter:on

    Matrix<double, 5, 6> abs_A = OmR * abs_L * OmL;

    QState _qL(zL); // Консервативный вектор слева
    QState _qR(zR); // Консервативный вектор справа

#ifdef ZEPHYR_ENABLE_EIGEN
    Vector6d qL, qR;
    qL.segment(0, 5) = _qL.vec();
    qL[5] = zL.pressure;
    qR.segment(0, 5) = _qR.vec();
    qR[5] = zR.pressure;

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа

    Flux res = 0.5 * (fL.vec() + fR.vec()) + 0.5 * abs_A * (qL - qR);

    return res;
#else
    throw std::runtime_error("CIR2 error: can't compute without eigen");
#endif
}

} // namespace zephyr