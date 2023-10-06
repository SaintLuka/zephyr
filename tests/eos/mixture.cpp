#include <iostream>
#include <iomanip>
#include <algorithm>

#include <zephyr/utils/error_list.h>

#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/eos/stiffened_gas.h>
#include <zephyr/phys/eos/mie_gruneisen.h>

#include <zephyr/phys/eos/materials.h>

using namespace zephyr::phys;
using namespace zephyr::utils;

/// @brief Проверяет производные (dE/dR)_P и (dE/dP)_R
/// (dV/dP)_T и (dV/dT)_P
/// (dE/dP)_T и (dE/dT)_P
/// Также проверяется аппроксимация двучленным УрС и скорость звука
void test_mixture(Materials& mixture, Fractions& beta, double rho, double eps) {
    std::cout << std::scientific << std::setprecision(2);

    dRdE P = mixture.pressure_re(rho, eps, beta, {.deriv = true});
    auto T = mixture.temperature_rp(rho, P, beta);
    dPdT V = mixture.volume_pt(P, T, beta);
    dPdT E = mixture.energy_pt(P, T, beta);

    double delta = 1.0e-5;
    double dR = delta * rho;
    double dE = delta * eps;
    double dT = delta * T;
    double dP = delta * P;

    double dPdR_E = (mixture.pressure_re(rho + dR, eps, beta) - mixture.pressure_re(rho - dR, eps, beta)) / (2 * dR);
    double dPdE_R = (mixture.pressure_re(rho, eps + dE, beta) - mixture.pressure_re(rho, eps - dE, beta)) / (2 * dE);
    double dVdP_T = (mixture.volume_pt(P + dP, T, beta) - mixture.volume_pt(P - dP, T, beta)) / (2 * dP);
    double dVdT_P = (mixture.volume_pt(P, T + dT, beta) - mixture.volume_pt(P, T - dT, beta)) / (2 * dT);
    double dEdP_T = (mixture.energy_pt(P + dP, T, beta) - mixture.energy_pt(P - dP, T, beta)) / (2 * dP);
    double dEdT_P = (mixture.energy_pt(P, T + dT, beta) - mixture.energy_pt(P, T - dT, beta)) / (2 * dT);

    ErrorList err1 = {
            {eps,  E.val},    // Совпадают величины
            {rho,  1.0 / V.val},

            {P.dR, dPdR_E},  // Совпадают производные
            {P.dE, dPdE_R},
            {V.dP, dVdP_T},
            {V.dT, dVdT_P},
            {E.dP, dEdP_T},
            {E.dT, dEdT_P}
    };

    if (err1.is_ok(1.0e-9)) {
        std::cout << "\tDerivatives.  Max error: " << err1.max() << ": OK\n";
    }
    else {
        std::cerr << "\tDerivatives.  Max error: " << err1.max() << "\n";
        throw std::runtime_error("Mixture test failed #1");
    }

    // Проверим аппроксимацию двучленным УрС
    StiffenedGas sg = mixture.stiffened_gas(rho, P, beta);
    dRdE P2 = sg.pressure_re(rho, eps, {.deriv= true});

    ErrorList err2 = {
            {P.val, P2.val},
            {P.dR,  P2.dR},
            {P.dE,  P2.dE}
    };

    if (err2.is_ok(1.0e-14)) {
        std::cout << "\tStiffenedGas. Max error: " << err2.max() << ": OK\n";
    } else {
        std::cerr << "\tStiffenedGas. Max error: " << err2.max() << "\n";
        throw std::runtime_error("Mixture test failed #2");
    }

    double c1 = std::sqrt(P.dR + P.val * P.dE / (rho * rho));
    double c2 = mixture.sound_speed_re(rho, eps, beta);
    double c3 = mixture.sound_speed_rp(rho, P.val, beta);
    double c4 = sg.sound_speed_re(rho, eps);

    ErrorList err3 = {
            {c1, c2},
            {c1, c3},
            {c1, c4}
    };

    if (err3.is_ok(1.0e-14)) {
        std::cout << "\tSound Speed.  Max error: " << err3.max() << ": OK\n";
    }
    else {
        std::cerr << "\tSound Speed.  Max error: " << err3.max() << "\n";
        throw std::runtime_error("Mixture test failed #3");
    }
}

struct State {
    ///< Максимальное число материалов
    static const int max_size = Fractions::max_size;

    double P;    ///< Равновесное давление
    double T;    ///< Равновесная температура
    double rho;  ///< Смесевая плотность
    double eps;  ///< Смесевая энергия

    Fractions alpha;  ///< Объемные концентрации
    Fractions beta;   ///< Массовые концентрации
    std::array<double, max_size> rs;  ///< Истинные плотности
    std::array<double, max_size> es;  ///< Истинные энергии

    /// @brief Вывести информацию о смеси
    /// @param n Актуальное число материалов
    void print(int n) const {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "\tTemperature: " << std::setw(8) << T - 0.0_C    << "  C\n";
        std::cout << "\tPressure:    " << std::setw(8) << 1.0e-6 * P   << "  MPa\n";
        std::cout << "\tMix density: " << std::setw(8) << 1.0e-3 * rho << "  g/cm^3\n";
        std::cout << "\tMix energy:  " << std::setw(8) << 1.0e-6 * eps << "  MJ/kg\n\n";

        std::cout << "\tvol. frac: ";
        for (int i = 0; i < n; ++i) {
            std::cout << std::setw(10) << 100.0 * alpha[i];
        }
        std::cout << "  %\n";

        std::cout << "\tmass frac: ";
        for (int i = 0; i < n; ++i) {
            std::cout << std::setw(10) << 100.0 * beta[i];
        }
        std::cout << "  %\n";

        std::cout << "\tdensity:   ";
        for (int i = 0; i < n; ++i) {
            std::cout << std::setw(10) << 1.0e-3 * rs[i];
        }
        std::cout << "  g/cm^3\n";

        std::cout << "\tenergy:    ";
        for (int i = 0; i < n; ++i) {
            std::cout << std::setw(10) << 1.0e-6 * es[i];
        }
        std::cout << "  MJ/kg\n";
    }
};

/// @brief Явные функции смеси
/// @param mixture Смесь материалов
/// @param alpha Объемные доли
/// @param P Давление смеси
/// @param T Температура смеси
void test_explicit(Materials& mixture, Fractions& alpha, double P, double T) {
    State info;

    info.P = P;
    info.T = T;
    info.alpha = alpha;

    int n = mixture.size();

    for (int i = 0; i < n; ++i) {
        info.rs[i] = 1.0 / mixture[i].volume_pt(P, T);
        info.es[i] = mixture[i].energy_pt(P, T);
    }

    info.rho = 0.0;
    for (int i = 0; i < n; ++i) {
        info.rho += alpha[i] * info.rs[i];
    }

    for (int i = 0; i < n; ++i) {
        info.beta[i] = alpha[i] * info.rs[i] / info.rho;
    }

    info.eps = mixture.energy_pt(P, T, info.beta);

    info.print(n);
}

/// @brief Неявные функции смеси
/// @param mixture Смесь материалов
/// @param beta Массовые доли
/// @param rho Смесевая плотность
/// @param P Равновесное давление
void test_implicit_rp(Materials& mixture, Fractions& beta, double rho, double P) {
    State info;

    info.P = P;
    info.rho = rho;
    info.beta = beta;

    int n = mixture.size();

    info.T = mixture.temperature_rp(rho, P, beta);
    info.eps = mixture.energy_pt(P, info.T, beta);

    for (int i = 0; i < n; ++i) {
        info.rs[i] = 1.0 / mixture[i].volume_pt(P, info.T);
        info.es[i] = mixture[i].energy_pt(P, info.T);

        info.alpha[i] = beta[i] * rho / info.rs[i];
    }

    info.print(n);
}

/// @brief Неявные функции смеси
/// @param mixture Смесь материалов
/// @param beta Массовые доли
/// @param rho Смесевая плотность
/// @param p Равновесное давление
void test_implicit_rt(Materials& mixture, Fractions& beta, double rho, double T) {
    State info;

    info.T = T;
    info.rho = rho;
    info.beta = beta;

    int n = mixture.size();

    info.P = mixture.pressure_rt(rho, T, beta);
    info.eps = mixture.energy_pt(info.P, info.T, beta);

    for (int i = 0; i < n; ++i) {
        info.rs[i] = 1.0 / mixture[i].volume_pt(info.P, info.T);
        info.es[i] = mixture[i].energy_pt(info.P, info.T);

        info.alpha[i] = beta[i] * rho / info.rs[i];
    }

    info.print(n);
}

/// @brief Неявные функции смеси
/// @param mixture Смесь материалов
/// @param beta Массовые доли
/// @param rho Смесевая плотность
/// @param eps Энергия смеси
void test_implicit_re(Materials& mixture, Fractions& beta, double rho, double eps) {
    State info;

    info.eps = eps;
    info.rho = rho;
    info.beta = beta;

    int n = mixture.size();

    info.P = mixture.pressure_re(rho, eps, beta);
    info.T = mixture.temperature_rp(rho, info.P, beta);

    for (int i = 0; i < n; ++i) {
        info.rs[i] = 1.0 / mixture[i].volume_pt(info.P, info.T);
        info.es[i] = mixture[i].energy_pt(info.P, info.T);

        info.alpha[i] = beta[i] * rho / info.rs[i];
    }

    info.print(n);
}

int main() {
    // Смесь
    Materials mixture;
    mixture += IdealGas::create("Air");
    mixture += StiffenedGas::create("Water2");
    mixture += MieGruneisen::create("Cu");

    std::cout << "Mixture:\n";
    std::cout << "\tMaterial[0]: IdealGas(\"Air\")\n";
    std::cout << "\tMaterial[1]: StiffenedGas(\"Water2\")\n";
    std::cout << "\tMaterial[2]: MieGruneisen(\"Cu\")\n";

    // Параметры (не согласованы, и не обязаны)
    double rho_test = 2.13_g_cm3;
    double eps_test = 0.29_MJ_kg;
    double P_test = 90.0_MPa;
    double T_test = 400.0_C;

    // Объемные доли
    Fractions alpha = {0.5, 0.3, 0.2};

    // Массовые концентрации
    Fractions beta = {0.11, 0.07, 0.82};


    std::cout << "\nTest mixture:\n";
    test_mixture(mixture, beta, rho_test, eps_test);

    std::cout << "\nTest mixture (explicit)\n";
    test_explicit(mixture, alpha, P_test, T_test);

    std::cout << "\nTest mixture (implicit, (rho, p))\n";
    test_implicit_rp(mixture, beta, rho_test, P_test);

    std::cout << "\nTest mixture (implicit, (rho, T))\n";
    test_implicit_rt(mixture, beta, rho_test, T_test);

    std::cout << "\nTest mixture (implicit, (rho, e))\n";
    test_implicit_re(mixture, beta, rho_test, eps_test);

    return 0;
}