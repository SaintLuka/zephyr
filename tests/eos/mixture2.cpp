#include <iostream>
#include <iomanip>
#include <algorithm>

#include <zephyr/geom/vector.h>

#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/eos/mie_gruneisen.h>
#include <zephyr/phys/eos/materials.h>

#include <zephyr/utils/matplotlib.h>

using namespace zephyr::geom;
using namespace zephyr::phys;
using namespace zephyr::utils;

namespace plt = zephyr::utils::matplotlib;

double fraction(double x) {
    double delta = 0.6;
    if (x > +delta) {
        return 0.0;
    } else if (x < -delta) {
        return 1.0;
    } else {
        double c = std::sin(0.25 * M_PI * (x + delta) / delta);
        return 1.0 - c * c;
    }
}

int main() {
    // Смесь
    auto eos1 = MieGruneisen::create("Pb");
    auto eos2 = IdealGas::create("Air");

    Materials mixture;
    mixture += eos1;
    mixture += eos2;

    auto x = linspace(-1.0, 1.0, 10000);
    //auto x = linspace(-0.65, -0.55, 1000);
    //auto x = linspace(-0.6, -0.595, 10000);

    std::vector<double> a1(x.size());
    std::vector<double> a2(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        a1[i] = fraction(x[i]);
        a2[i] = 1.0 - a1[i];
    }

    // Истинные значения
    double dens1 = 10100.0;
    double pres1 = -1.0e9;
    double vol1  = 1.0 / dens1;
    double enrg1 = eos1->energy_rp(dens1, pres1);
    double temp1 = eos1->temperature_rp(dens1, pres1);
    double dens2 = 5.0;
    double pres2 = 8.0e6;
    double vol2  = 1.0 / dens2;
    double enrg2 = eos2->energy_rp(dens2, pres2);
    double tempA = eos2->temperature_rp(dens2, pres2);

    std::vector<double> mix_dens(x.size());
    std::vector<double> mix_energy(x.size());
    std::vector<double> mix_pres(x.size());
    std::vector<double> mix_temp(x.size());
    std::vector<double> beta1(x.size());
    std::vector<double> beta2(x.size());

    std::vector<double> mix_gamma(x.size());
    std::vector<double> mix_p_min(x.size());
    std::vector<double> mix_e0(x.size());

    for (size_t i = 0; i < x.size(); ++i) {
        mix_dens[i] = a1[i] * dens1 + a2[i] * dens2;

        beta1[i] = a1[i] * dens1 / mix_dens[i];
        beta2[i] = a2[i] * dens2 / mix_dens[i];

        mix_energy[i] = beta1[i] * enrg1 + beta2[i] * enrg2;

        Fractions beta = {beta1[i], beta2[i]};
        mix_pres[i] = mixture.pressure_re(mix_dens[i], mix_energy[i], beta);
        mix_temp[i] = mixture.temperature_rp(mix_dens[i], mix_pres[i], beta);

        StiffenedGas gas = mixture.stiffened_gas(mix_dens[i], mix_pres[i], beta);
        mix_gamma[i] = gas.gamma;
        mix_p_min[i] = gas.min_pressure();
        mix_e0[i] = gas.eps_0;
    }

    plt::figure_size(15.0, 8.0);

    plt::subplot2grid(2, 3, 0, 0);
    plt::title("Vol Fractions");
    //plt::plot(x, a1, "k-o");
    plt::plot(x, a2);

    plt::subplot2grid(2, 3, 1, 0);
    plt::title("Mass Fractions");
    //plt::plot(x, beta1);
    plt::plot(x, beta2);

    plt::subplot2grid(2, 3, 0, 1);
    plt::title("Mix density");
    plt::plot(x, mix_dens);

    plt::subplot2grid(2, 3, 0, 2);
    plt::title("Mix energy");
    plt::plot(x, mix_energy);

    plt::subplot2grid(2, 3, 1, 1);
    plt::title("Mix pressure");
    plt::plot(x, mix_pres);

    plt::subplot2grid(2, 3, 1, 2);
    plt::title("Mix temperature");
    plt::plot(x, mix_temp);

    plt::tight_layout();
    plt::show();

    return 0;
}