#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <zephyr/utils/error_list.h>

#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/eos/stiffened_gas.h>
#include <zephyr/phys/eos/mie_gruneisen.h>

using namespace zephyr::phys;
using namespace zephyr::utils;


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
            {T,   eos.temperature_rp(rho, P1)}
    };

    if (err1.is_ok(1.0e-14)) {
        std::cout << "\tConsistency.  Max error: " << err1.max() << ": OK!\n";
    } else {
        std::cerr << "\tConsistency.  Max error: " << err1.max() << "\n";
        throw std::runtime_error("Eos test failed #1");
    }

    // Выражения с производными
    auto P = eos.pressure_re(rho, eps, {.deriv = true});
    auto V = eos.volume_pt(P, T, {.deriv = true});
    auto E = eos.energy_pt(P, T, {.deriv = true});

    // Проверка скорости звука
    double c = std::sqrt(P.dR + P.val * P.dE / (rho * rho));

    ErrorList err2 = {
            {c, eos.sound_speed_rp(rho, P)},
            {c, eos.sound_speed_re(rho, E)}
    };

    if (err2.is_ok(1.0e-14)) {
        std::cout << "\tSound Speed.  Max error: " << err2.max() << ": OK!\n";
    } else {
        std::cerr << "\tSound Speed.  Max error: " << err2.max() << "\n";
        throw std::runtime_error("Eos test failed #2");
    }

    // Аппроксимация двучленным УрС
    StiffenedGas sg = eos.stiffened_gas(rho, P);

    ErrorList err3 = {
            {c,    sg.sound_speed_re(rho, E)},
            {c,    sg.sound_speed_rp(rho, P)},
            {P,    sg.pressure_re(rho, E)},
            {P.dR, sg.pressure_re(rho, E, {.deriv = true}).dR},
            {P.dE, sg.pressure_re(rho, E, {.deriv = true}).dE}
    };

    if (err3.is_ok(1.0e-14)) {
        std::cout << "\tStiffenedGas. Max error: " << err3.max() << ": OK!\n";
    } else {
        std::cerr << "\tStiffenedGas. Max error: " << err3.max() << "\n";
        throw std::runtime_error("Eos test failed #3");
    }

    // Проверка производных
    double d = 1.0e-5;

    double dP = d * P;
    double dR = d * rho;
    double dT = d * T;
    double dE = d * eps;

    double dPdR_e = (eos.pressure_re(rho + dR, eps) - eos.pressure_re(rho - dR, eps)) / (2 * dR);
    double dPdE_r = (eos.pressure_re(rho, eps + dE) - eos.pressure_re(rho, eps - dE)) / (2 * dE);
    double dVdP_t = (eos.volume_pt(P + dP, T) - eos.volume_pt(P - dP, T)) / (2 * dP);
    double dVdT_p = (eos.volume_pt(P, T + dT) - eos.volume_pt(P, T - dT)) / (2 * dT);
    double dEdP_t = (eos.energy_pt(P + dP, T) - eos.energy_pt(P - dP, T)) / (2 * dP);
    double dEdT_p = (eos.energy_pt(P, T + dT) - eos.energy_pt(P, T - dT)) / (2 * dT);

    ErrorList err4 = {
            {P.dR, dPdR_e},
            {P.dE, dPdE_r},
            {V.dP, dVdP_t},
            {V.dT, dVdT_p},
            {E.dP, dEdP_t},
            {E.dT, dEdT_p},
    };

    if (err4.is_ok(1.0e-9)) {
        std::cout << "\tDerivatives.  Max error: " << err4.max() << ": OK!\n";
    } else {
        std::cerr << "\tDerivatives.  Max error: " << err4.max() << "\n";
        throw std::runtime_error("Eos test failed #4");
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