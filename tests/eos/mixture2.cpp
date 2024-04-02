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
    auto air = IdealGas::create("Air");
    auto lead = MieGruneisen::create("Pb");

    Materials mixture;
    mixture += lead;
    mixture += air;

    auto x = linspace(-1.0, 1.0, 10000);

    std::vector<double> a1(x.size());
    std::vector<double> a2(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        a1[i] = fraction(x[i]);
        a2[i] = 1.0 - a1[i];
    }

    // Истинные значения
    double densL = 10100.0;
    double presL = -1.0e9;
    double tempL = lead->temperature_rp(densL, presL);
    double densA = 5.0;
    double presA = 8.0e6;
    double tempA = air->temperature_rp(densA, presA);

    std::vector<double> mix_dens(x.size());
    std::vector<double> mix_energy(x.size());
    std::vector<double> mix_pres(x.size());
    std::vector<double> mix_temp(x.size());

    for (size_t i = 0; i < x.size(); ++i) {
        mix_dens[i] = 1.0 / (a1[i] * lead->volume_pt(presL, tempL) + a2[i] * air->volume_pt(presA, tempA));
        mix_energy[i] = a1[i] * lead->energy_pt(presL, tempL) + a2[i] * air->energy_pt(presA, tempA);

        Fractions beta = {a1[i] * densL / mix_dens[i], a2[i] * densA / mix_dens[i] };
        mix_pres[i] = mixture.pressure_re(mix_dens[i], mix_energy[i], beta);
        mix_temp[i] = mixture.temperature_rp(mix_dens[i], mix_pres[i], beta);
    }

    plt::figure_size(10.0, 8.0);

    plt::subplot2grid(2, 2, 0, 0);
    plt::title("Mix density");
    plt::plot(x, mix_dens);

    plt::subplot2grid(2, 2, 0, 1);
    plt::title("Mix energy");
    plt::plot(x, mix_energy);

    plt::subplot2grid(2, 2, 1, 0);
    plt::title("Mix pressure");
    plt::plot(x, mix_pres);

    plt::subplot2grid(2, 2, 1, 1);
    plt::title("Mix temperature");
    plt::plot(x, mix_temp);

    plt::tight_layout();
    plt::show();

    return 0;
}