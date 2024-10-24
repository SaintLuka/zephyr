#include <iostream>
#include <iomanip>
#include <algorithm>

#include <zephyr/utils/error_list.h>

#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/eos/stiffened_gas.h>
#include <zephyr/phys/eos/mie_gruneisen.h>

#include <zephyr/phys/eos/materials.h>

#include <zephyr/utils/matplotlib.h>
namespace plt = zephyr::utils::matplotlib;

using namespace zephyr::phys;
using namespace zephyr::utils;


int main() {
    // Смесь
    auto Cu = MieGruneisen::create("Cu");
    auto Air = IdealGas::create("Air");
    Materials mixture;
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

    auto rho_tests = zephyr::geom::linspace(0.9*rho_test, 1.1*rho_test, 1000);
    std::vector<double> Pcu_rb0(rho_tests.size());
    std::vector<double> Pcu_rb1(rho_tests.size());
    std::vector<double> Pcu_rb2(rho_tests.size());
    std::vector<double> Pcu_rb3(rho_tests.size());

    std::vector<double> Pcu_rb0_r(rho_tests.size());
    std::vector<double> Pcu_rb1_r(rho_tests.size());
    std::vector<double> Pcu_rb2_r(rho_tests.size());
    std::vector<double> Pcu_rb3_r(rho_tests.size());

    std::vector<double> Pcu_rb0_e(rho_tests.size());
    std::vector<double> Pcu_rb1_e(rho_tests.size());
    std::vector<double> Pcu_rb2_e(rho_tests.size());
    std::vector<double> Pcu_rb3_e(rho_tests.size());

    for (size_t i = 0; i < rho_tests.size(); ++i) {
        Pcu_rb0[i] = mixture.pressure_re(rho_tests[i], e_test, {1.0, 0.0});
        Pcu_rb1[i] = mixture.pressure_re(rho_tests[i], e_test, {1.0 - sb1, sb1});
        Pcu_rb2[i] = mixture.pressure_re(rho_tests[i], e_test, {1.0 - sb2, sb2});
        Pcu_rb3[i] = mixture.pressure_re(rho_tests[i], e_test, {1.0 - sb3, sb3});

        auto [P0, T0] = mixture.find_PT(rho_tests[i], e_test, {1.0, 0.0});

        Pcu_rb0_r[i] = mixture.pressure_re(rho_tests[i], e_test, {1.0, 0.0}, {.deriv=true, .P0=P0, .T0=T0}).dR;
        Pcu_rb1_r[i] = mixture.pressure_re(rho_tests[i], e_test, {1.0 - sb1, sb1}, {.deriv=true, .P0=P0, .T0=T0}).dR;
        Pcu_rb2_r[i] = mixture.pressure_re(rho_tests[i], e_test, {1.0 - sb2, sb2}, {.deriv=true, .P0=P0, .T0=T0}).dR;
        Pcu_rb3_r[i] = mixture.pressure_re(rho_tests[i], e_test, {1.0 - sb3, sb3}, {.deriv=true, .P0=P0, .T0=T0}).dR;

        Pcu_rb0_e[i] = mixture.pressure_re(rho_tests[i], e_test, {1.0, 0.0}, {.deriv=true, .P0=P0, .T0=T0}).dE;
        Pcu_rb1_e[i] = mixture.pressure_re(rho_tests[i], e_test, {1.0 - sb1, sb1}, {.deriv=true, .P0=P0, .T0=T0}).dE;
        Pcu_rb2_e[i] = mixture.pressure_re(rho_tests[i], e_test, {1.0 - sb2, sb2}, {.deriv=true, .P0=P0, .T0=T0}).dE;
        Pcu_rb3_e[i] = mixture.pressure_re(rho_tests[i], e_test, {1.0 - sb3, sb3}, {.deriv=true, .P0=P0, .T0=T0}).dE;
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



    auto e_tests   = zephyr::geom::linspace(0.95*e_test, 1.05*e_test, 100);
    std::vector<double> Pcu_eb0(e_tests.size());
    std::vector<double> Pcu_eb1(e_tests.size());
    std::vector<double> Pcu_eb2(e_tests.size());
    std::vector<double> Pcu_eb3(e_tests.size());

    std::vector<double> Pcu_eb0_r(e_tests.size());
    std::vector<double> Pcu_eb1_r(e_tests.size());
    std::vector<double> Pcu_eb2_r(e_tests.size());
    std::vector<double> Pcu_eb3_r(e_tests.size());

    std::vector<double> Pcu_eb0_e(e_tests.size());
    std::vector<double> Pcu_eb1_e(e_tests.size());
    std::vector<double> Pcu_eb2_e(e_tests.size());
    std::vector<double> Pcu_eb3_e(e_tests.size());

    std::vector<double> Pcu_eb0_c(e_tests.size());
    std::vector<double> Pcu_eb1_c(e_tests.size());
    std::vector<double> Pcu_eb2_c(e_tests.size());
    std::vector<double> Pcu_eb3_c(e_tests.size());

    for (size_t i = 0; i < e_tests.size(); ++i) {
        Pcu_eb0[i] = mixture.pressure_re(rho_test, e_tests[i], {1.0, 0.0});
        Pcu_eb1[i] = mixture.pressure_re(rho_test, e_tests[i], {1.0 - sb1, sb1});
        Pcu_eb2[i] = mixture.pressure_re(rho_test, e_tests[i], {1.0 - sb2, sb2});
        Pcu_eb3[i] = mixture.pressure_re(rho_test, e_tests[i], {1.0 - sb3, sb3});

        auto [P0, T0] = mixture.find_PT(rho_test, e_tests[i], {1.0, 0.0});

        Pcu_eb0_r[i] = mixture.pressure_re(rho_test, e_tests[i], {1.0, 0.0}, {.deriv=true, .P0=P0, .T0=T0}).dR;
        Pcu_eb1_r[i] = mixture.pressure_re(rho_test, e_tests[i], {1.0 - sb1, sb1}, {.deriv=true, .P0=P0, .T0=T0}).dR;
        Pcu_eb2_r[i] = mixture.pressure_re(rho_test, e_tests[i], {1.0 - sb2, sb2}, {.deriv=true, .P0=P0, .T0=T0}).dR;
        Pcu_eb3_r[i] = mixture.pressure_re(rho_test, e_tests[i], {1.0 - sb3, sb3}, {.deriv=true, .P0=P0, .T0=T0}).dR;

        Pcu_eb0_e[i] = mixture.pressure_re(rho_test, e_tests[i], {1.0, 0.0}, {.deriv=true, .P0=P0, .T0=T0}).dE;
        Pcu_eb1_e[i] = mixture.pressure_re(rho_test, e_tests[i], {1.0 - sb1, sb1}, {.deriv=true, .P0=P0, .T0=T0}).dE;
        Pcu_eb2_e[i] = mixture.pressure_re(rho_test, e_tests[i], {1.0 - sb2, sb2}, {.deriv=true, .P0=P0, .T0=T0}).dE;
        Pcu_eb3_e[i] = mixture.pressure_re(rho_test, e_tests[i], {1.0 - sb3, sb3}, {.deriv=true, .P0=P0, .T0=T0}).dE;
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