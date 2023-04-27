#include <iostream>

#include <zephyr/geom/vector.h>
#include <zephyr/math/cfd/flux/hllc.h>

namespace zephyr {
namespace math {

using geom::Vector5d;
using geom::Vector6d;
using geom::Matrix5d;


smf::Flux HLLC2::flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux HLLC2::calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) {
    using namespace smf;

    double rho1 = zL.density, rho2 = zR.density;
    double u1 = zL.velocity.x(), u2 = zR.velocity.x();
    double v1 = zL.velocity.y(), v2 = zR.velocity.y();
    double w1 = zL.velocity.z(), w2 = zR.velocity.z();
    double p1 = zL.pressure, p2 = zR.pressure;

    double c1 = eos.sound_speed_rp(rho1, p1); // скорость звука слева std::sqrt(gamma * pressure / density);
    double c2 = eos.sound_speed_rp(rho2, p2); // скорость звука справа

    double s1 = std::min(u1 - c1, u2 - c2); // SL
    double s2 = std::max(u1 + c1, u2 + c2); // SR
    double S = (p2 - p1 + rho1 * u1 * (s1 - u1) - rho2 * u2 * (s2 - u2)) / (rho1 * (s1 - u1) - rho2 * (s2 - u2)); // S*

    double energy_state1 = zL.energy / rho1 + (S - u1) * (S + p1 / rho1 / (s1 - u1));
    double energy_state2 = zR.energy / rho2 + (S - u2) * (S + p2 / rho2 / (s2 - u2));

    double coeff1 = rho1 * (s1 - u1) / (s1 - S);
    double coeff2 = rho2 * (s2 - u2) / (s2 - S);

    QState QL(coeff1, Vector3d(coeff1 * S, coeff1 * v1, coeff1 * w1), coeff1 * energy_state1);
    QState QR(coeff2, Vector3d(coeff2 * S, coeff2 * v2, coeff2 * w2), coeff2 * energy_state2);

    QState qL(zL); // Консервативный вектор слева
    QState qR(zR); // Консервативный вектор справа

    Flux fL(zL);   // Дифференциальный поток слева
    Flux fR(zR);   // Дифференциальный поток справа

    Flux F_hllc;
    if (0 <= s1)
        F_hllc = fL;
    else if (s1 <= 0 && 0 <= S)
        F_hllc = fL.vec() + s1 * (QL.vec() - qL.vec());
    else if (S <= 0 && 0 <= s2)
        F_hllc = fR.vec() + s2 * (QR.vec() - qR.vec());
    else if (s2 <= 0)
        F_hllc = fR;
    else {
        std::cerr << "zL: " << zL << "\n";
        std::cerr << "Zr: " << zR << "\n";
        std::cerr << "Sound speed left: " << c1 << " , Sound speed right: " << c2 << "\n";
        std::cerr << "SL: " << s1 << ", SR: " << s2 << ", S*: " << S << "\n";
        throw std::runtime_error("HLLC:calc_flux Error, strange case in switch");
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

    double c1 = eos.sound_speed_rp(rho1, p1); // скорость звука слева std::sqrt(gamma * pressure / density);
    double c2 = eos.sound_speed_rp(rho2, p2); // скорость звука справа

    double s1 = std::min(u1 - c1, u2 - c2); // SL
    double s2 = std::max(u1 + c1, u2 + c2); // SR
    double S = (p2 - p1 + rho1 * u1 * (s1 - u1) - rho2 * u2 * (s2 - u2)) / (rho1 * (s1 - u1) - rho2 * (s2 - u2)); // S*

    double rhoL = rho1 * (s1 - u1) / (s1 - S); // rho*_L
    double rhoR = rho2 * (s2 - u2) / (s2 - S); // rho*_R

//    double P = p1 + rho1 * (s1 - u1) / (S - u1); // P* (на 0 часто делится, пока не знаю что с этим делать)
    double P = (p1 + p2) / 2; // P* альтернативная формула
//    double P = (p1 * rho2 * c2 + p2 * rho1 * c1 + (u1 - u2) * rho1 * rho2 * c1 * c2) / (rho1 * c1 + rho2 * c2); ещё одна альтернативная формула

    double VL = (s1 - u1) / (s1 - S) * v1; // V*_L
    double VR = (s2 - u2) / (s2 - S) * v2; // V*_R
    double WL = (s1 - u1) / (s1 - S) * w1; // W*_L
    double WR = (s2 - u2) / (s2 - S) * w2; // W*_R
    double EL = eos.energy_rp(rhoL, P);
    double ER = eos.energy_rp(rhoR, P);

    PState stateL(rhoL, Vector3d(S, VL, WL), P, EL);
    PState stateR(rhoR, Vector3d(S, VR, WR), P, ER);

    QState qL(zL); // U_L
    QState qR(zR); // U_R
    QState QL(stateL); // U*_L
    QState QR(stateR); // U*_R

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
        throw std::runtime_error("HLLC2::calc_flux Error, strange case in switch");
    }

    if (std::isnan(F_hllc.mass) || std::isnan(F_hllc.momentum.x()) || std::isnan(F_hllc.momentum.y()) ||
        std::isnan(F_hllc.momentum.z()) || std::isnan(F_hllc.energy)) {
        std::cerr << "zL: " << zL << "\n";
        std::cerr << "Zr: " << zR << "\n";
        std::cerr << "Sound speed left: " << c1 << " , Sound speed right: " << c2 << "\n";
        std::cerr << "SL: " << s1 << ", SR: " << s2 << ", S*: " << S << "\n";
        std::cerr << "F_HLLLC: " << F_hllc << "\n";
        throw std::runtime_error("HLLC2::calc_flux Error, F_HLLC has the bad value");
    }

    return F_hllc;
}

}
}