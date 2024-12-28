/// @file Проверка совместности для уравнений состояния

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <zephyr/math/calc/derivatives.h>
#include <zephyr/utils/error_list.h>

#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/matter/eos/stiffened_gas.h>
#include <zephyr/phys/matter/eos/mie_gruneisen.h>

#include <zephyr/utils/matplotlib.h>
namespace plt = zephyr::utils::matplotlib;


using namespace zephyr::phys;
using namespace zephyr::utils;
using namespace zephyr::math;


void test_eos(Eos& eos, double T, double rho) {
    double P1 = eos.pressure_rT(rho, T);
    double eps = eos.energy_rP(rho, P1);

    // Проверка совместности
    ErrorList err1 = {
            {P1,  eos.pressure_re(rho, eps)},
            {P1,  eos.pressure_re(rho, eps)},
            {eps, eos.energy_rP(rho, P1)},
            {P1,  eos.pressure_rT(rho, T)},
            {rho, 1.0 / eos.volume_PT(P1, T)},
            {eps, eos.energy_PT(P1, T)},
            {eps, eos.energy_rT(P1, T)},
            {T,   eos.temperature_rP(rho, P1)}
    };

    if (err1.is_ok(1.0e-14)) {
        std::cout << "\tConsistency.  Max error: " << err1.max() << ": OK!\n";
    } else {
        std::cout << "\tConsistency.  Max error: " << err1.max() << ": NOT OK!!\n";
    }

    // Выражения с производными
    auto P = eos.pressure_re(rho, eps, {.deriv = true});
    auto V = eos.volume_PT(P, T, {.deriv = true});
    auto E1 = eos.energy_PT(P, T, {.deriv = true});
    auto E2 = eos.energy_rT(rho, T, {.deriv = true});


    // Проверка скорости звука
    double c = std::sqrt(P.dR + P.val * P.dE / (rho * rho));

    ErrorList err2 = {
            {c, eos.sound_speed_rP(rho, P)},
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
            {c,    sg.sound_speed_rP(rho, P)},
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
                return eos.volume_PT(x, T);
            }, P, dP);

    double dVdT_p = derivative<1, 4>(
            [&](double x) -> double {
                return eos.volume_PT(P, x);
            }, T, dT);

    double dEdP_t = derivative<1, 4>(
            [&](double x) -> double {
                return eos.energy_PT(x, T);
            }, P, dP);

    double dEdT_p = derivative<1, 4>(
            [&](double x) -> double {
                return eos.energy_PT(P, x);
            }, T, dT);

    double dEdR_t = derivative<1, 4>(
            [&](double x) -> double {
                return eos.energy_rT(x, T);
            }, rho, dR);

    double dEdT_r = derivative<1, 4>(
            [&](double x) -> double {
                return eos.energy_rT(rho, x);
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
    std::cout << std::setprecision(16);

    StiffenedGas sg("HeavyGas");

    std::cout << sg.temperature_rP(4.0, 1.0) << "\n";

    auto rho = zephyr::geom::linspace(100, 15000, 200);

    for (auto T: {0.0, 300.0, 800.0}) {
        std::vector<double> P1(rho.size());
        for (int i = 0; i < rho.size(); ++i) {
            P1[i] = sg.pressure_rT(rho[i], T);
        }

        plt::plot(rho, P1, {{"color", "blue"}, {"label", "T" + std::to_string(T)}});
    }

    plt::show();



    /*
    std::cout << "IdealGas(\"Air\")\n";
    IdealGas eos1("Air");
    test_eos(eos1, 145.0_C, 1.3_kg_m3);

    std::cout << "\nStiffenedGas(\"Water2\")\n";
    StiffenedGas eos2("Water2");
    test_eos(eos2, 65.0_C, 1.1_g_cm3);

    std::cout << "\nMieGruneisen(\"Cu\")\n";
    MieGruneisen eos3("Cu");
    test_eos(eos3, 470.0_C, 9.1_g_cm3);
     */

    return 0;
}