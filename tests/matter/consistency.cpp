// @brief Проверка совместности для уравнений состояния

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <zephyr/math/calc/derivatives.h>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/matter/eos/stiffened_gas.h>
#include <zephyr/phys/matter/eos/mie_gruneisen.h>

#include <zephyr/utils/error_list.h>

using namespace zephyr::math;
using namespace zephyr::phys;
using namespace zephyr::utils;

struct State {
    double rho, e, T, P;

    State(const Eos& eos, double density, double temperature) {
        rho = density;
        T = temperature;
        e = eos.energy_rT(rho, T);
        P = eos.pressure_rT(rho, T);
    }

    void print() const {
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "\tTemperature: " << std::setw(10) << (T - 0.0_C)    << "  C\n";
        std::cout << "\tPressure:    " << std::setw(10) << (1.0e-6 * P) << "  MPa\n";
        std::cout << "\tDensity:     " << std::setw(10) << (1.0e-3 * rho) << "  g/cm^3\n";
        std::cout << "\tEnergy:      " << std::setw(10) << (1.0e-6 * e)   << "  MJ/kg\n\n";
    }
};

// Проверка совместности
void test_consistency(const Eos& eos, State z) {
    std::cout << std::scientific << std::setprecision(2);

    double P1 = eos.pressure_re(z.rho, z.e);
    double P2 = eos.pressure_rT(z.rho, z.T);
    double e1 = eos.energy_rP(z.rho, z.P);
    double e2 = eos.energy_rT(z.rho, z.T);
    double T1 = eos.temperature_rP(z.rho, z.P);
    double e3 = eos.energy_PT(z.P, z.T);
    double v1 = eos.volume_PT(z.P, z.T);

    ErrorList err = {
            {z.P,   P1},
            {z.P,   P2},
            {z.e,   e1},
            {z.e,   e2},
            {z.e,   e3},
            {z.T,   T1},
            {z.rho, 1.0 / v1}
    };

    if (err.is_ok(1.0e-14)) {
        std::cout << "\tConsistency.  Max error: " << std::setw(12) << err.max() << ".\tOK!\n";
    } else {
        std::cout << "\tConsistency.  Max error: " << std::setw(12) << err.max() << ".\tNOT OK!!\n";
        err.print("\t\t");
    }
}

// Выражения с производными
void test_derivatives(const Eos& eos, State z) {
    double d = 1.0e-3;
    double dP = d * z.P;
    double dR = d * z.rho;
    double dT = d * z.T;
    double dE = d * z.e;

    ErrorList err;

    // P(rho, e)
    {
        auto P = eos.pressure_re(z.rho, z.e, {.deriv = true});

        double dPdR_E = derivative<1, 4>([&](double rho) -> double {
            return eos.pressure_re(rho, z.e);
        }, z.rho, dR);
        double dPdE_R = derivative<1, 4>([&](double e) -> double {
            return eos.pressure_re(z.rho, e);
        }, z.e, dE);

        err += {P.dR, dPdR_E};
        err += {P.dE, dPdE_R};
    }

    // P(rho, T)
    {
        auto P = eos.pressure_rT(z.rho, z.T, {.deriv = true});

        double dPdR_T = derivative<1, 4>([&](double rho) -> double {
            return eos.pressure_rT(rho, z.T);
        }, z.rho, dR);
        double dPdT_R = derivative<1, 4>([&](double T) -> double {
            return eos.pressure_rT(z.rho, T);
        }, z.T, dT);

        err += {P.dR, dPdR_T};
        err += {P.dT, dPdT_R};
    }

    // e(rho, P)
    {
        auto E = eos.energy_rP(z.rho, z.P, {.deriv = true});

        double dEdR_P = derivative<1, 4>([&](double rho) -> double {
            return eos.energy_rP(rho, z.P);
        }, z.rho, dR);
        double dEdP_R = derivative<1, 4>([&](double P) -> double {
            return eos.energy_rP(z.rho, P);
        }, z.P, dP);

        err += {E.dR, dEdR_P};
        err += {E.dP, dEdP_R};
    }

    // e(rho, T)
    {
        auto E = eos.energy_rT(z.rho, z.T, {.deriv = true});

        double dEdR_T = derivative<1, 4>([&](double rho) -> double {
            return eos.energy_rT(rho, z.T);
        }, z.rho, dR);
        double dEdT_R = derivative<1, 4>([&](double T) -> double {
            return eos.energy_rT(z.rho, T);
        }, z.T, dT);

        err += {E.dR, dEdR_T};
        err += {E.dT, dEdT_R};
    }

    // v(P, T)
    {
        auto V = eos.volume_PT(z.P, z.T, {.deriv = true});
        double dVdP_T = derivative<1, 4>([&](double P) -> double {
            return eos.volume_PT(P, z.T);
        }, z.P, dP);
        double dVdT_P = derivative<1, 4>([&](double T) -> double {
            return eos.volume_PT(z.P, T);
        }, z.T, dT);

        err += {V.dP, dVdP_T};
        err += {V.dT, dVdT_P};
    }

    // e(P, T)
    {
        auto E = eos.energy_PT(z.P, z.T, {.deriv = true});
        double dEdP_T = derivative<1, 4>([&](double P) -> double {
            return eos.energy_PT(P, z.T);
        }, z.P, dP);
        double dEdT_P = derivative<1, 4>([&](double T) -> double {
            return eos.energy_PT(z.P, T);
        }, z.T, dT);

        err += {E.dP, dEdP_T};
        err += {E.dT, dEdT_P};
    }

    if (err.is_ok(1.0e-10)) {
        std::cout << "\tDerivatives.  Max error: " << std::setw(12) << err.max() << ".\tOK!\n";
    } else {
        std::cout << "\tDerivatives.  Max error: " << std::setw(12) << err.max() << ".\tNOT OK!!\n";
        err.print("\t\t");
    }
}

// Определение скорости звука
void test_sound_speed(const Eos& eos, State z) {
    dRdE P = eos.pressure_re(z.rho, z.e, {.deriv=true});
    dRdP e = eos.energy_rP  (z.rho, z.P, {.deriv=true});

    // Проверка скорости звука
    double c1 = std::sqrt(P.dR + P.val * P.dE / (z.rho * z.rho));
    double c2 = std::sqrt((P / (z.rho * z.rho) - e.dR) / e.dP);
    double c3 = eos.sound_speed_rP(z.rho, z.P);
    double c4 = eos.sound_speed_re(z.rho, z.e);

    ErrorList err = {
            {c1, c2},
            {c1, c3},
            {c1, c4}
    };

    if (err.is_ok(1.0e-14)) {
        std::cout << "\tSound Speed.  Max error: " << std::setw(12) << err.max() << ".\tOK!\n";
    } else {
        std::cout << "\tSound Speed.  Max error: " << std::setw(12) << err.max() << ".\tNOT OK!!\n";
        err.print("\t\t");
    }
}

// Аппроксимация двучленным УрС
void test_stiffened_gas(const Eos& eos, State z) {
    StiffenedGas sg = eos.stiffened_gas(z.rho, z.P);

    dRdE P = eos.pressure_re(z.rho, z.e, {.deriv=true});
    dRdP e = eos.energy_rP  (z.rho, z.P, {.deriv=true});
    double c = eos.sound_speed_rP(z.rho, z.P);

    ErrorList err = {
            {c,    sg.sound_speed_re(z.rho, z.e)},
            {c,    sg.sound_speed_rP(z.rho, z.P)},
            {P,    sg.pressure_re(z.rho, z.e)},
            {P.dR, sg.pressure_re(z.rho, z.e, {.deriv = true}).dR},
            {P.dE, sg.pressure_re(z.rho, z.e, {.deriv = true}).dE},
            {e,    sg.energy_rP(z.rho, z.P)},
            {e.dR, sg.energy_rP(z.rho, z.P, {.deriv = true}).dR},
            {e.dP, sg.energy_rP(z.rho, z.P, {.deriv = true}).dP}
    };

    if (err.is_ok(1.0e-13)) {
        std::cout << "\tStiffenedGas. Max error: " << std::setw(12) << err.max() << ".\tOK!\n";
    } else {
        std::cout << "\tStiffenedGas. Max error: " << std::setw(12) << err.max() << ".\tNOT OK!!\n";
        err.print("\t\t");
    }
}

void test_eos(Eos& eos, double density, double temperature) {
    State z(eos, density, temperature);
    z.print();

    test_consistency(eos, z);
    test_derivatives(eos, z);
    test_sound_speed(eos, z);
    test_stiffened_gas(eos, z);
}

int main() {
    std::cout << std::setprecision(16);
    std::cout << "IdealGas(\"Air\")\n";
    IdealGas eos1("Air");
    test_eos(eos1, 1.3_kg_m3, 145.0_C);

    std::cout << "\nStiffenedGas(\"Water2\")\n";
    StiffenedGas eos2("Water2");
    test_eos(eos2, 1.1_g_cm3, 65.0_C);

    std::cout << "\nMieGruneisen(\"Cu\")\n";
    MieGruneisen eos3("Cu");
    test_eos(eos3, 9.1_g_cm3, 470.0_C);

    return 0;
}