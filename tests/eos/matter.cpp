#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <zephyr/math/calc/derivatives.h>
#include <zephyr/utils/error_list.h>

#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/eos/strange_gas.h>
#include <zephyr/phys/eos/stiffened_gas.h>
#include <zephyr/phys/eos/mie_gruneisen.h>


using namespace zephyr::phys;
using namespace zephyr::utils;
using namespace zephyr::math;


void test_eos(Eos& eos, double T, double rho) {
    double P1 = eos.pressure_rt(rho, T);
    double eps = eos.energy_rp(rho, P1);

    // Проверка совместности
    ErrorList err1 = {
            {P1,  eos.pressure_re(rho, eps)},
            {P1,  eos.pressure_re(rho, eps)},
            {eps, eos.energy_rp(rho, P1)},
            {P1,  eos.pressure_rt(rho, T)},
            {rho, 1.0 / eos.volume_pt(P1, T)},
            {eps, eos.energy_pt(P1, T)},
            {eps, eos.energy_rt(P1, T)},
            {T,   eos.temperature_rp(rho, P1)}
    };

    if (err1.is_ok(1.0e-14)) {
        std::cout << "\tConsistency.  Max error: " << err1.max() << ": OK!\n";
    } else {
        std::cout << "\tConsistency.  Max error: " << err1.max() << ": NOT OK!!\n";
    }

    // Выражения с производными
    auto P = eos.pressure_re(rho, eps, {.deriv = true});
    auto V = eos.volume_pt(P, T, {.deriv = true});
    auto E1 = eos.energy_pt(P, T, {.deriv = true});
    auto E2 = eos.energy_rt(rho, T, {.deriv = true});


    // Проверка скорости звука
    double c = std::sqrt(P.dR + P.val * P.dE / (rho * rho));

    ErrorList err2 = {
            {c, eos.sound_speed_rp(rho, P)},
            {c, eos.sound_speed_re(rho, E1)}
    };

    if (err2.is_ok(1.0e-14)) {
        std::cout << "\tSound Speed.  Max error: " << err2.max() << ": OK!\n";
    } else {
        std::cout << "\tSound Speed.  Max error: " << err2.max() << ": NOT OK!!\n";
    }

    // Аппроксимация двучленным УрС
    StiffenedGas sg = eos.stiffened_gas(rho, P);

    ErrorList err3 = {
            {c,    sg.sound_speed_re(rho, E1)},
            {c,    sg.sound_speed_rp(rho, P)},
            {P,    sg.pressure_re(rho, E1)},
            {P.dR, sg.pressure_re(rho, E1, {.deriv = true}).dR},
            {P.dE, sg.pressure_re(rho, E1, {.deriv = true}).dE}
    };

    if (err3.is_ok(1.0e-13)) {
        std::cout << "\tStiffenedGas. Max error: " << err3.max() << ": OK!\n";
    } else {
        std::cout << "\tStiffenedGas. Max error: " << err3.max() << ": NOT OK!!\n";
        std::cout << err3 << "\n";
    }

    // Проверка производных
    double d = 1.0e-3;

    double dP = d * P;
    double dR = d * rho;
    double dT = d * T;
    double dE = d * eps;

    double dPdR_e = derivative<1, 4>(
            [&](double x) -> double {
                return eos.pressure_re(x, eps);
            }, rho, dR);

    double dPdE_r = derivative<1, 4>(
            [&](double x) -> double {
                return eos.pressure_re(rho, x);
            }, eps, dE);

    double dVdP_t = derivative<1, 4>(
            [&](double x) -> double {
                return eos.volume_pt(x, T);
            }, P, dP);

    double dVdT_p = derivative<1, 4>(
            [&](double x) -> double {
                return eos.volume_pt(P, x);
            }, T, dT);

    double dEdP_t = derivative<1, 4>(
            [&](double x) -> double {
                return eos.energy_pt(x, T);
            }, P, dP);

    double dEdT_p = derivative<1, 4>(
            [&](double x) -> double {
                return eos.energy_pt(P, x);
            }, T, dT);

    double dEdR_t = derivative<1, 4>(
            [&](double x) -> double {
                return eos.energy_rt(x, T);
            }, rho, dR);

    double dEdT_r = derivative<1, 4>(
            [&](double x) -> double {
                return eos.energy_rt(rho, x);
            }, T, dT);

    ErrorList err4 = {
            {P.dR,  dPdR_e},
            {P.dE,  dPdE_r},
            {V.dP,  dVdP_t},
            {V.dT,  dVdT_p},
            {E1.dP, dEdP_t},
            {E1.dT, dEdT_p},
            {E2.dR, dEdR_t},
            {E2.dT, dEdT_r},
    };

    if (err4.is_ok(1.0e-8)) {
        std::cout << "\tDerivatives.  Max error: " << err4.max() << ": OK!\n";
    } else {
        std::cout << "\tDerivatives.  Max error: " << err4.max() << ": NOT OK!!\n";
        std::cout << err4 << "\n";
    }
}

int main() {
    std::cout << std::scientific << std::setprecision(2);

    std::cout << "IdealGas(\"Air\")\n";
    IdealGas eos1("Air");
    test_eos(eos1, 145.0_C, 1.3_kg_m3);

    std::cout << "\nStiffenedGas(\"Water2\")\n";
    StiffenedGas eos2("Water2");
    test_eos(eos2, 65.0_C, 1.1_g_cm3);

    std::cout << "\nMieGruneisen(\"Cu\")\n";
    MieGruneisen eos3("Cu");
    test_eos(eos3, 470.0_C, 9.1_g_cm3);

    return 0;
}