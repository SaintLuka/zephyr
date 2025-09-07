#include <iostream>
#include <iomanip>
#include <algorithm>

#include <zephyr/utils/error_list.h>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/matter/eos/stiffened_gas.h>
#include <zephyr/phys/matter/eos/mie_gruneisen.h>

#include <zephyr/phys/matter/mixture_pt.h>

#include <zephyr/utils/numpy.h>
#include <zephyr/utils/range.h>
#include <zephyr/utils/matplotlib.h>

namespace plt = zephyr::utils::matplotlib;

using namespace zephyr;
using namespace zephyr::phys;
using namespace zephyr::utils;


int main() {
    // Смесь
    auto Cu = MieGruneisen::create("Cu");
    auto Air = IdealGas::create("Air");
    MixturePT mixture;
    mixture += Cu;
    mixture += Air;

    std::cout << "Mixture:\n";
    std::cout << "\tMaterial[0]: MieGruneisen(\"Cu\")\n";
    std::cout << "\tMaterial[1]: IdealGas(\"Air\")\n\n";

    // Параметры (не согласованы, и не обязаны)
    double rho_test = 9.1_g_cm3;
    double P_test = 1.0e9;
    double T_test = Cu->temperature_rP(rho_test, P_test);
    double e_test = Cu->energy_rP(rho_test, P_test);

    std::cout << "dens_mix: " << rho_test << "\n";
    std::cout << "energy_mix: " << e_test << "\n";
    std::cout << "P_test: " << P_test << "\n";
    std::cout << "T_test: " << T_test << "\n\n";

    double sb1 = 1.0e-8;
    double sb2 = 1.0e-6;
    double sb3 = 1.0e-4;

    plt::figure_size(16.0, 8.0, 150);

    auto rho_tests = np::linspace(0.9*rho_test, 1.1*rho_test, 1000);
    auto Pcu_rb0 = np::empty_like(rho_tests);
    auto Pcu_rb1 = np::empty_like(rho_tests);
    auto Pcu_rb2 = np::empty_like(rho_tests);
    auto Pcu_rb3 = np::empty_like(rho_tests);

    auto Pcu_rb0_r = np::empty_like(rho_tests);
    auto Pcu_rb1_r = np::empty_like(rho_tests);
    auto Pcu_rb2_r = np::empty_like(rho_tests);
    auto Pcu_rb3_r = np::empty_like(rho_tests);

    auto Pcu_rb0_e = np::empty_like(rho_tests);
    auto Pcu_rb1_e = np::empty_like(rho_tests);
    auto Pcu_rb2_e = np::empty_like(rho_tests);
    auto Pcu_rb3_e = np::empty_like(rho_tests);

    for (auto [i, rho]: enumerate(rho_tests)) {
        Pcu_rb0[i] = mixture.pressure_re(rho, e_test, {1.0, 0.0});
        Pcu_rb1[i] = mixture.pressure_re(rho, e_test, {1.0 - sb1, sb1});
        Pcu_rb2[i] = mixture.pressure_re(rho, e_test, {1.0 - sb2, sb2});
        Pcu_rb3[i] = mixture.pressure_re(rho, e_test, {1.0 - sb3, sb3});

        auto [rhos, P0, T0] = mixture.get_rPT(rho, e_test, {1.0, 0.0});

        Pcu_rb0_r[i] = mixture.pressure_re(rho, e_test, {1.0, 0.0}, {.deriv=true, .P0=P0, .T0=T0}).dR;
        Pcu_rb1_r[i] = mixture.pressure_re(rho, e_test, {1.0 - sb1, sb1}, {.deriv=true, .P0=P0, .T0=T0}).dR;
        Pcu_rb2_r[i] = mixture.pressure_re(rho, e_test, {1.0 - sb2, sb2}, {.deriv=true, .P0=P0, .T0=T0}).dR;
        Pcu_rb3_r[i] = mixture.pressure_re(rho, e_test, {1.0 - sb3, sb3}, {.deriv=true, .P0=P0, .T0=T0}).dR;

        Pcu_rb0_e[i] = mixture.pressure_re(rho, e_test, {1.0, 0.0}, {.deriv=true, .P0=P0, .T0=T0}).dE;
        Pcu_rb1_e[i] = mixture.pressure_re(rho, e_test, {1.0 - sb1, sb1}, {.deriv=true, .P0=P0, .T0=T0}).dE;
        Pcu_rb2_e[i] = mixture.pressure_re(rho, e_test, {1.0 - sb2, sb2}, {.deriv=true, .P0=P0, .T0=T0}).dE;
        Pcu_rb3_e[i] = mixture.pressure_re(rho, e_test, {1.0 - sb3, sb3}, {.deriv=true, .P0=P0, .T0=T0}).dE;
    }

    plt::subplot2grid(2, 3, 0, 0);
    plt::xlabel("Density");
    plt::ylabel("Pressure");
    plt::plot(rho_tests, Pcu_rb0, "g");
    plt::plot(rho_tests, Pcu_rb1, "k--");
    plt::plot(rho_tests, Pcu_rb2, "k--");
    plt::plot(rho_tests, Pcu_rb3, "k--");

    plt::subplot2grid(2, 3, 0, 1);
    plt::xlabel("Density");
    plt::ylabel("$dP/dR \\, ( \\approx c^2 )$");
    plt::plot(rho_tests, Pcu_rb0_r, "g");
    plt::plot(rho_tests, Pcu_rb1_r, "k--");
    plt::plot(rho_tests, Pcu_rb2_r, "k--");
    plt::plot(rho_tests, Pcu_rb3_r, "k--");

    plt::subplot2grid(2, 3, 0, 2);
    plt::xlabel("Density");
    plt::ylabel("dP/dE");
    plt::plot(rho_tests, Pcu_rb0_e, "g");
    plt::plot(rho_tests, Pcu_rb1_e, "k--");
    plt::plot(rho_tests, Pcu_rb2_e, "k--");
    plt::plot(rho_tests, Pcu_rb3_e, "k--");


    auto e_tests = np::linspace(0.95*e_test, 1.05*e_test, 100);
    auto Pcu_eb0 = np::empty_like(e_tests);
    auto Pcu_eb1 = np::empty_like(e_tests);
    auto Pcu_eb2 = np::empty_like(e_tests);
    auto Pcu_eb3 = np::empty_like(e_tests);

    auto Pcu_eb0_r = np::empty_like(e_tests);
    auto Pcu_eb1_r = np::empty_like(e_tests);
    auto Pcu_eb2_r = np::empty_like(e_tests);
    auto Pcu_eb3_r = np::empty_like(e_tests);

    auto Pcu_eb0_e = np::empty_like(e_tests);
    auto Pcu_eb1_e = np::empty_like(e_tests);
    auto Pcu_eb2_e = np::empty_like(e_tests);
    auto Pcu_eb3_e = np::empty_like(e_tests);

    auto Pcu_eb0_c = np::empty_like(e_tests);
    auto Pcu_eb1_c = np::empty_like(e_tests);
    auto Pcu_eb2_c = np::empty_like(e_tests);
    auto Pcu_eb3_c = np::empty_like(e_tests);

    for (auto [i, engy]: enumerate(e_tests)) {
        Pcu_eb0[i] = mixture.pressure_re(rho_test, engy, {1.0, 0.0});
        Pcu_eb1[i] = mixture.pressure_re(rho_test, engy, {1.0 - sb1, sb1});
        Pcu_eb2[i] = mixture.pressure_re(rho_test, engy, {1.0 - sb2, sb2});
        Pcu_eb3[i] = mixture.pressure_re(rho_test, engy, {1.0 - sb3, sb3});

        auto [rhos, P0, T0] = mixture.get_rPT(rho_test, engy, {1.0, 0.0});

        Pcu_eb0_r[i] = mixture.pressure_re(rho_test, engy, {1.0, 0.0}, {.deriv=true, .P0=P0, .T0=T0}).dR;
        Pcu_eb1_r[i] = mixture.pressure_re(rho_test, engy, {1.0 - sb1, sb1}, {.deriv=true, .P0=P0, .T0=T0}).dR;
        Pcu_eb2_r[i] = mixture.pressure_re(rho_test, engy, {1.0 - sb2, sb2}, {.deriv=true, .P0=P0, .T0=T0}).dR;
        Pcu_eb3_r[i] = mixture.pressure_re(rho_test, engy, {1.0 - sb3, sb3}, {.deriv=true, .P0=P0, .T0=T0}).dR;

        Pcu_eb0_e[i] = mixture.pressure_re(rho_test, engy, {1.0, 0.0}, {.deriv=true, .P0=P0, .T0=T0}).dE;
        Pcu_eb1_e[i] = mixture.pressure_re(rho_test, engy, {1.0 - sb1, sb1}, {.deriv=true, .P0=P0, .T0=T0}).dE;
        Pcu_eb2_e[i] = mixture.pressure_re(rho_test, engy, {1.0 - sb2, sb2}, {.deriv=true, .P0=P0, .T0=T0}).dE;
        Pcu_eb3_e[i] = mixture.pressure_re(rho_test, engy, {1.0 - sb3, sb3}, {.deriv=true, .P0=P0, .T0=T0}).dE;
    }

    plt::subplot2grid(2, 3, 1, 0);
    plt::xlabel("Energy");
    plt::ylabel("Pressure");
    plt::plot(e_tests, Pcu_eb0, "g");
    plt::plot(e_tests, Pcu_eb1, "k--");
    plt::plot(e_tests, Pcu_eb2, "k--");
    plt::plot(e_tests, Pcu_eb3, "k--");

    plt::subplot2grid(2, 3, 1, 1);
    plt::xlabel("Energy");
    plt::ylabel("$dP/dR \\, ( \\approx c^2 )$");
    plt::plot(e_tests, Pcu_eb0_r, "g");
    plt::plot(e_tests, Pcu_eb1_r, "k--");
    plt::plot(e_tests, Pcu_eb2_r, "k--");
    plt::plot(e_tests, Pcu_eb3_r, "k--");

    plt::subplot2grid(2, 3, 1, 2);
    plt::xlabel("Energy");
    plt::ylabel("dP/dE");
    plt::plot(e_tests, Pcu_eb0_e, "g");
    plt::plot(e_tests, Pcu_eb1_e, "k--");
    plt::plot(e_tests, Pcu_eb2_e, "k--");
    plt::plot(e_tests, Pcu_eb3_e, "k--");


    plt::tight_layout();
    plt::show();

    return 0;
}