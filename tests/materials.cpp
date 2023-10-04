#include <iostream>
#include <iomanip>

#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/eos/stiffened_gas.h>

#include <zephyr/phys/eos/materials.h>

using namespace zephyr::phys;

struct State {
    ///< Максимальное число материалов
    static const int max_size = Fractions::max_size;

    double p;    ///< Равновесное давление
    double T;    ///< Равновесная температура
    double rho;  ///< Смесевая плотность
    double e;    ///< Смесевая энергия

    Fractions alpha;  ///< Объемные концентрации
    Fractions beta;   ///< Массовые концентрации
    std::array<double, max_size> rhos;  ///< Истинные плотности
    std::array<double, max_size> es;    ///< Истинные энергии

    /// @brief Вывести информацию о смеси
    /// @param n Актуальное число материалов
    void print(int n) const {
        //std::cout << std::fixed << std::setprecision(1) << "\n";
        for (int i = 0; i < n; ++i) {
            std::cout << "Material(" << i << ")\n";
            std::cout << "  vol frac:  " << 100.0 * alpha[i] << "%\n";
            std::cout << "  mass frac: " << 100.0 * beta[i] << "%\n";
            std::cout << "  density:   " << rhos[i] << "\n";
            std::cout << "  energy:    " << es[i] << "\n";
        }

        std::cout << "Pressure:    " << p << "\n";
        std::cout << "Temperature: " << T << "\n";
        std::cout << "Mix density: " << rho << "\n";
        std::cout << "Mix energy:  " << e << "\n";
    }
};

/// @brief Явные функции смеси
/// @param mixture Смесь материалов
/// @param alpha Объемные доли
/// @param p Давление смеси
/// @param T Температура смеси
void test_explicit(Materials& mixture, Fractions& alpha, double p, double T) {
    State info;

    info.p = p;
    info.T = T;
    info.alpha = alpha;

    int n = mixture.size();

    for (int i = 0; i < n; ++i) {
        info.rhos[i] = 1.0 / mixture[i].volume_pt(p, T);
        info.es[i] = mixture[i].energy_pt(p, T);
    }

    info.rho = 0.0;
    for (int i = 0; i < n; ++i) {
        info.rho += alpha[i] * info.rhos[i];
    }

    for (int i = 0; i < n; ++i) {
        info.beta[i] = alpha[i] * info.rhos[i] / info.rho;
    }

    info.e = mixture.energy_pt(p, T, info.beta);

    info.print(n);
}

/// @brief Неявные функции смеси
/// @param mixture Смесь материалов
/// @param beta Массовые доли
/// @param rho Смесевая плотность
/// @param p Равновесное давление
void test_implicit_rp(Materials& mixture, Fractions& beta, double rho, double p) {
    State info;

    info.p = p;
    info.rho = rho;
    info.beta = beta;

    int n = mixture.size();

    info.T = mixture.temperature_rp(rho, p, beta);
    info.e = mixture.energy_pt(p, info.T, beta);

    for (int i = 0; i < n; ++i) {
        info.rhos[i] = 1.0 / mixture[i].volume_pt(p, info.T);
        info.es[i] = mixture[i].energy_pt(p, info.T);

        info.alpha[i] = beta[i] * rho / info.rhos[i];
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

    info.p = mixture.pressure_rt(rho, T, beta);
    info.e = mixture.energy_pt(info.p, info.T, beta);

    for (int i = 0; i < n; ++i) {
        info.rhos[i] = 1.0 / mixture[i].volume_pt(info.p, info.T);
        info.es[i] = mixture[i].energy_pt(info.p, info.T);

        info.alpha[i] = beta[i] * rho / info.rhos[i];
    }

    info.print(n);
}

/// @brief Неявные функции смеси
/// @param mixture Смесь материалов
/// @param beta Массовые доли
/// @param rho Смесевая плотность
/// @param e Энергия смеси
void test_implicit_re(Materials& mixture, Fractions& beta, double rho, double e) {
    State info;

    info.e = e;
    info.rho = rho;
    info.beta = beta;

    int n = mixture.size();

    info.p = mixture.pressure_re(rho, e, beta);
    info.T = mixture.temperature_rp(rho, info.p, beta);

    for (int i = 0; i < n; ++i) {
        info.rhos[i] = 1.0 / mixture[i].volume_pt(info.p, info.T);
        info.es[i] = mixture[i].energy_pt(info.p, info.T);

        info.alpha[i] = beta[i] * rho / info.rhos[i];
    }

    info.print(n);
}

int main() {
    std::cout << "Test eos\n";

    Eos::test();

    std::cout << "OK!\n\n";


    std::cout << "Test mixture (explicit)\n";

    // Смесь
    Materials mixture;
    mixture += IdealGas::create("Air");
    mixture += StiffenedGas::create("Water2");

    // Объемные доли
    Fractions alpha = {0.9, 0.1};

    test_explicit(mixture, alpha, 1.0_bar, 20.0_C);

    // Массовые концентрации
    Fractions beta = {0.01, 0.99};

    std::cout << "\nTest mixture (implicit, (rho, p))\n";
    test_implicit_rp(mixture, beta, 105.0, 1.0_bar);

    std::cout << "\nTest mixture (implicit, (rho, T))\n";
    test_implicit_rt(mixture, beta, 210.0, 20.0_C);

    std::cout << "\nTest mixture (implicit, (rho, e))\n";
    test_implicit_re(mixture, beta, 115.0, 90.0e3);

    return 0;
}