#include <string>
#include <iostream>
#include <stdexcept>
#include <cmath>

#include <zephyr/phys/eos/mie_gruneisen.h>
#include <zephyr/phys/eos/stiffened_gas.h>

namespace zephyr { namespace phys {

MieGruneisen::MieGruneisen(const std::string &name) {
    table_params(name);
}

double MieGruneisen::cold_pressure(double rho) const {
    return B / n * (std::pow(v0 * rho, n) - 1.0);
}

double MieGruneisen::cold_energy(double rho) const {
    return B / (n * (n - 1.0)) * ((std::pow(v0 * rho, n) - 1.0 + n) / rho - n * v0);
}

dRdE MieGruneisen::pressure_re(double rho, double eps, const Options& options) const {
    double P_c = cold_pressure(rho);
    double e_c = cold_energy(rho);
    double P = P_c + Gr * rho * (eps - e_c);
    dRdE res{P};
    if (options.deriv) {
        res.dR = (P + B + (n - Gr - 1.0) * P_c) / rho;
        res.dE = Gr * rho;
    }
    return res;
}

double MieGruneisen::energy_rp(double rho, double P, const Options& options) const {
    double P_c = cold_pressure(rho);
    double e_c = cold_energy(rho);
    return e_c + (P - P_c) / (Gr * rho);
}

double MieGruneisen::sound_speed_re(double rho, double eps, const Options& options) const {
    double P_c = cold_pressure(rho);
    double e_c = cold_energy(rho);
    double c2 = Gr * (1.0 + Gr) * (eps - e_c) + (B + n * P_c) / rho;
    return std::sqrt(std::max(c2, 1.0e-6));
}

double MieGruneisen::sound_speed_rp(double rho, double P, const Options& options) const {
    double P_c = cold_pressure(rho);
    double c2 = (B + (1.0 + Gr) * P + (n - Gr - 1.0) * P_c) / rho;
    return std::sqrt(std::max(c2, 1.0e-6));
}

double MieGruneisen::pressure_rt(double rho, double T, const Options& options) const {
    double P_c = cold_pressure(rho);
    return P_c + Gr * Cv * rho * (T - T0);
}

double MieGruneisen::temperature_rp(double rho, double P, const Options& options) const {
    double P_c = cold_pressure(rho);
    return T0 + (P - P_c) / (Gr * Cv * rho);
}

// Решить уравнение A * x - x^(-nu) = C, где nu > 0
// При A > 0 существует единственное решение,
// условие нарушается при P < p_min = - B/n
double solve(double A, double C, double nu, double x0) {
    double x = x0;
    double err = 1.0;
    int counter = 0;
    while (err > 1.0e-12 && counter < 30) {
        double pow = std::pow(x, -nu);
        double dx = (C - A * x + pow) / (A + nu * pow / x);
        err = std::abs(dx / x);
        x += dx;
        ++counter;
    }
    return x;
}

dPdT MieGruneisen::volume_pt(double P, double T, const Options& options) const {
    if (P < min_pressure()) {
        return NAN;
    }

    double A = 1.0 + n * P / B;
    double nu = n - 1.0;
    double C = n * Gr * Cv * (T - T0) / (B * v0);

    double v_init = std::isnan(options.rho0) ? 1.0 : 1.0 / (v0 * options.rho0);

    double vol = v0 * solve(A, C, nu, v_init);
    dPdT res{vol};
    if (options.deriv) {
        double P_c = cold_pressure(1.0 / vol);
        double xi = 1.0 / (B + P + (n - 1.0) * P_c);
        res.dP = -vol * xi;
        res.dT = Gr * Cv * xi;
    }
    return res;
}

dPdT MieGruneisen::energy_pt(double P, double T, const Options& options) const {
    if (P < min_pressure()) {
        return NAN;
    }

    dPdT vol = volume_pt(P, T, options);
    double rho = 1.0 / vol;
    dPdT res{energy_rp(rho, P, options)};
    if (options.deriv) {
        double mP_c = -cold_pressure(rho);
        res.dP = mP_c * vol.dP;
        res.dT = mP_c * vol.dT + Cv;
    }
    return res;
}

StiffenedGas MieGruneisen::stiffened_gas(double rho, double P, const Options& options) const {
    double P_c = cold_pressure(rho);
    double e_c = cold_energy(rho);

    double gamma = 1.0 + Gr;
    double eps_0 = e_c + ((Gr - n) * P_c - B) / (Gr * rho);
    double P0 = (B + (n - Gr - 1.0) * P_c) / gamma;

    return StiffenedGas(gamma, P0, eps_0, NAN);
}

double MieGruneisen::min_pressure() const {
    return -B / n;
}

void MieGruneisen::table_params(std::string name) {
    double r0, c0;
    T0 = 0.0_C;

    // @formatter:off
    if (name == "Ag") {
        r0 = 10490.0;
        c0 = 3178.0;
        n  = 4.5;
        Gr = 2.55;
        Cv = 227.0;
    } else if (name == "Air") {
        r0 = 1.16;
        c0 = 0.001;
        n  = 2.0;
        Gr = 0.4;
        Cv = 718.390;
    } else if (name == "Al") {
        r0 = 2710.0;
        c0 = 5333.0;
        n  = 3.5;
        Gr = 2.13;
        Cv = 865.0;
    } else if (name == "Au") {
        r0 = 19300.0;
        c0 = 3063.0;
        n  = 3.4;
        Gr = 3.29;
        Cv = 124.0;
    } else if (name == "Be") {
        r0 = 1850.0;
        c0 = 7993.0;
        n  = 3.0;
        Gr = 1.17;
        Cv = 1630;
    } else if (name == "Bi") {
        r0 = 9790.0;
        c0 = 1850.0;
        n  = 4.9;
        Gr = 0.17;
        Cv = 123;
    } else if (name == "Cd") {
        r0 = 8640.0;
        c0 = 2434.0;
        n  = 4.5;
        Gr = 2.32;
        Cv = 218.0;
    } else if (name == "Co") {
        r0 = 8820.0;
        c0 = 4743.0;
        n  = 3.6;
        Gr = 1.97;
        Cv = 414.0;
    } else if (name == "Cr") {
        r0 = 7100.0;
        c0 = 5153.0;
        n  = 4.3;
        Gr = 1.19;
        Cv = 444.0;
    } else if (name == "Cu") {
        r0 = 8930.0;
        c0 = 3899.0;
        n  = 4.0;
        Gr = 2.0;
        Cv = 392.0;
    } else if (name == "Fe") {
        r0 = 7855.0;
        c0 = 3883.0;
        n  = 4.5;
        Gr = 1.6;
        Cv = 441.0;
    } else if (name == "In") {
        r0 = 7270.0;
        c0 = 2430.0;
        n  = 4.0;
        Gr = 2.238;
        Cv = 218.0;
    } else if (name == "Mg") {
        r0 = 1740.0;
        c0 = 4540.0;
        n  = 3.3;
        Gr = 1.46;
        Cv = 984.0;
    } else if (name == "Mo") {
        r0 = 10200.0;
        c0 = 5100.0;
        n  = 3.3;
        Gr = 1.52;
        Cv = 247.0;
    } else if (name == "Nb") {
        r0 = 8604.0;
        c0 = 4472.0;
        n  = 3.0;
        Gr = 1.679;
        Cv = 267.0;
    } else if (name == "Ni") {
        r0 = 8870.0;
        c0 = 4501.0;
        n  = 3.0;
        Gr = 1.83;
        Cv = 435.0;
    } else if (name == "Pb") {
        r0 = 11340.0;
        c0 = 1981.0;
        n  = 3.7;
        Gr = 2.77;
        Cv = 121.0;
        T0 = 137.078_C;
    } else if (name == "Pd") {
        r0 = 11950.0;
        c0 = 3955.0;
        n  = 4.7;
        Gr = 2.183;
        Cv = 239.0;
    } else if (name == "Pt") {
        r0 = 21420.0;
        c0 = 3605.0;
        n  = 4.0;
        Gr = 2.627;
        Cv = 130.0;
    } else if (name == "Rh") {
        r0 = 12420.0;
        c0 = 4775.0;
        n  = 3.8;
        Gr = 2.265;
        Cv = 240.0;
    } else if (name == "Sn") {
        r0 = 7280.0;
        c0 = 2437.0;
        n  = 4.1;
        Gr = 2.11;
        Cv = 203.0;
    } else if (name == "Ta") {
        r0 = 16450.0;
        c0 = 3431.0;
        n  = 3.0;
        Gr = 1.689;
        Cv = 139.0;
    } else if (name == "Th") {
        r0 = 11680.0;
        c0 = 2323.0;
        n  = 2.7;
        Gr = 1.26;
        Cv = 116.0;
    } else if (name == "Ti") {
        r0 = 4522.0;
        c0 = 5003.0;
        n  = 2.5;
        Gr = 1.09;
        Cv = 218.0;
    } else if (name == "Tl") {
        r0 = 11860.0;
        c0 = 1809.0;
        n  = 3.9;
        Gr = 2.25;
        Cv = 122.0;
    } else if (name == "V") {
        r0 = 6080.0;
        c0 = 5071.0;
        n  = 3.2;
        Gr = 1.29;
        Cv = 481.0;
    } else if (name == "W") {
        r0 = 19170.0;
        c0 = 4015.0;
        n  = 3.3;
        Gr = 1.54;
        Cv = 132.0;
    } else if (name == "Zn") {
        r0 = 7140.0;
        c0 = 3031.0;
        n  = 4.2;
        Gr = 2.45;
        Cv = 369.0;
    } else if (name == "Zr") {
        r0 = 6490.0;
        c0 = 3787.0;
        n  = 4.4;
        Gr = 0.771;
        Cv = 280.0;
    } else {
        throw std::runtime_error("Unknown Mie-Gruneisen material '" + name + "'");
    }
    // @formatter:on

    B = r0 * c0 * c0;
    v0 = 1.0 / r0;
}

}
}