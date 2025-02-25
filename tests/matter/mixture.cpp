#include <iostream>
#include <iomanip>
#include <algorithm>

#include <zephyr/math/calc/derivatives.h>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/matter/eos/stiffened_gas.h>
#include <zephyr/phys/matter/eos/mie_gruneisen.h>

#include <zephyr/phys/matter/mixture_pt.h>

#include <zephyr/utils/error_list.h>

using namespace zephyr::math;
using namespace zephyr::phys;
using namespace zephyr::utils;

const double gam_1 = 1.4;
const double gam_2 = 1.6;
const double gam_3 = 1.9;
const double Cv_1 = 1.0;
const double Cv_2 = 2.5;
const double Cv_3 = 3.5;


/// @brief Термодинамически согласованное состояние смеси
struct State {
    double P   = NAN;  ///< Равновесное давление
    double T   = NAN;  ///< Равновесная температура
    double rho = NAN;  ///< Смесевая плотность
    double e   = NAN;  ///< Смесевая энергия

    Fractions alpha = Fractions::Zero();  ///< Объемные концентрации
    Fractions beta  = Fractions::Zero();  ///< Массовые концентрации
    ScalarSet rs    = ScalarSet::NaN();   ///< Истинные плотности
    ScalarSet es    = ScalarSet::NaN();   ///< Истинные энергии

    /// @brief Создать сместь от (P, T)
    static State Explicit(MixturePT &mixture, Fractions alpha, double P, double T) {
        State z = {.P = P, .T = T, .alpha = alpha};

        for (int i = 0; i < mixture.size(); ++i) {
            z.rs[i] = 1.0 / mixture[i].volume_PT(P, T);
            z.es[i] = mixture[i].energy_PT(P, T);
        }

        z.rho = 0.0;
        for (int i = 0; i < mixture.size(); ++i) {
            z.rho += alpha[i] * z.rs[i];
        }

        for (int i = 0; i < mixture.size(); ++i) {
            z.beta[i] = alpha[i] * z.rs[i] / z.rho;
        }

        z.e = mixture.energy_PT(P, T, z.beta);
        return z;
    }

    /// @brief Создать сместь от (rho, P)
    static State Implicit_rP(MixturePT &mixture, Fractions &beta, double rho, double P) {
        State z = {.P = P, .rho = rho, .beta = beta};

        z.T = mixture.temperature_rP(rho, P, beta);
        z.e = mixture.energy_PT(P, z.T, beta);

        for (int i = 0; i < mixture.size(); ++i) {
            z.rs[i] = 1.0 / mixture[i].volume_PT(P, z.T);
            z.es[i] = mixture[i].energy_PT(P, z.T);
            z.alpha[i] = beta[i] * rho / z.rs[i];
        }
        return z;
    }

    /// @brief Создать сместь от (rho, T)
    static State Implicit_rT(MixturePT &mixture, Fractions &beta, double rho, double T) {
        State z = {.T = T, .rho = rho, .beta = beta};

        z.P = mixture.pressure_rT(rho, T, beta);
        z.e = mixture.energy_PT(z.P, z.T, beta);

        for (int i = 0; i < mixture.size(); ++i) {
            z.rs[i] = 1.0 / mixture[i].volume_PT(z.P, z.T);
            z.es[i] = mixture[i].energy_PT(z.P, z.T);
            z.alpha[i] = beta[i] * rho / z.rs[i];
        }
        return z;
    }

    /// @brief Создать сместь от (rho, e)
    static State Implicit_re(MixturePT &mixture, Fractions &beta, double rho, double e) {
        State z = {.rho = rho, .e = e, .beta = beta};

        z.P = mixture.pressure_re(rho, e, beta);
        z.T = mixture.temperature_rP(rho, z.P, beta);

        for (int i = 0; i < mixture.size(); ++i) {
            z.rs[i] = 1.0 / mixture[i].volume_PT(z.P, z.T);
            z.es[i] = mixture[i].energy_PT(z.P, z.T);
            z.alpha[i] = beta[i] * rho / z.rs[i];
        }
        return z;
    }

    inline static double format(double x) {
        std::cout << (x < 0.01 ? std::scientific : std::fixed);
        return x;
    }

    /// @brief Вывести информацию о смеси
    void print() const {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "\tTemperature: " << std::setw(10) << format(T - 0.0_C)    << "  C\n";
        std::cout << "\tPressure:    " << std::setw(10) << format(1.0e-6 * P) << "  MPa\n";
        std::cout << "\tMix density: " << std::setw(10) << format(1.0e-3 * rho) << "  g/cm^3\n";
        std::cout << "\tMix energy:  " << std::setw(10) << format(1.0e-6 * e)   << "  MJ/kg\n\n";

        std::cout << "\tvol. frac: ";
        for (int i = 0; i < beta.count(); ++i) {
            std::cout << std::setw(12) << format(100.0 * alpha[i]);
        }
        std::cout << "  %\n";

        std::cout << "\tmass frac: ";
        for (int i = 0; i < beta.count(); ++i) {
            std::cout << std::setw(12) << format(100.0 * beta[i]);
        }
        std::cout << "  %\n";

        std::cout << "\tdensity:   ";
        for (int i = 0; i < beta.count(); ++i) {
            std::cout << std::setw(12) << format(1.0e-3 * rs[i]);
        }
        std::cout << "  g/cm^3\n";

        std::cout << "\tenergy:    ";
        for (int i = 0; i < beta.count(); ++i) {
            std::cout << std::setw(12) << format(1.0e-6 * es[i]);
        }
        std::cout << "  MJ/kg\n";
    }
};

/// @brief Проверяет функции УрС
void test_consistency(MixturePT& mix, const State&z) {
    std::cout << std::scientific << std::setprecision(2);

    double P1 = mix.pressure_re(z.rho, z.e, z.beta);
    double P2 = mix.pressure_rT(z.rho, z.T, z.beta);
    double e1 = mix.energy_rT(z.rho, z.T, z.beta);
    double e2 = mix.energy_rP(z.rho, z.P, z.beta);
    double T1 = mix.temperature_rP(z.rho, z.P, z.beta);
    double r1 = 1.0 / mix.volume_PT(z.P, z.T, z.beta);
    double e3 = mix.energy_PT(z.P, z.T, z.beta);

    auto[rhos1, P3, T2] = mix.get_rPT(r1, e1, z.beta);
    auto[rhos2, e4, T3] = mix.get_reT(r1, P1, z.beta);

    // Совпадают величины
    ErrorList err = {
            {z.rho, r1},
            {z.e,   e1},
            {z.e,   e2},
            {z.e,   e3},
            {z.e,   e4},
            {z.P,   P1},
            {z.P,   P2},
            {z.P,   P3},
            {z.T,   T1},
            {z.T,   T2},
            {z.T,   T3}
    };
    for (int i = 0; i < z.alpha.count(); ++i) {
        err += {z.rs[i], rhos1[i]};
        err += {z.rs[i], rhos2[i]};
    }

    if (err.is_ok(1.0e-12)) {
        std::cout << "\tConsistency.  Max error: " << std::setw(12) << err.max() << ".\tOK!\n";
    }
    else {
        std::cout << "\tConsistency.  Max error: " << std::setw(12) << err.max() << ".\tNOT OK!\n";
        err.print("\t\t");
    }
}

/// @brief Проверяет производные. Сравнение с численными
void test_derivatives(MixturePT& mix, State z) {
    std::cout << std::scientific << std::setprecision(2);

    double delta = 1.0e-3;
    double dR = delta * z.rho;
    double dE = delta * z.e;
    double dT = delta * z.T;
    double dP = delta * z.P;

    ErrorList err;

    // P(rho, e)
    {
        auto P = mix.pressure_re(z.rho, z.e, z.beta, {.deriv = true});

        double dPdR_E = derivative<1, 4>([&](double rho) -> double {
            return mix.pressure_re(rho, z.e, z.beta);
        }, z.rho, dR);
        double dPdE_R = derivative<1, 4>([&](double e) -> double {
            return mix.pressure_re(z.rho, e, z.beta);
        }, z.e, dE);

        err += {P.dR, dPdR_E};
        err += {P.dE, dPdE_R};
    }

    // P(rho, T)
    {
        auto P = mix.pressure_rT(z.rho, z.T, z.beta, {.deriv = true});

        double dPdR_T = derivative<1, 4>([&](double rho) -> double {
            return mix.pressure_rT(rho, z.T, z.beta);
        }, z.rho, dR);
        double dPdT_R = derivative<1, 4>([&](double T) -> double {
            return mix.pressure_rT(z.rho, T, z.beta);
        }, z.T, dT);

        err += {P.dR, dPdR_T};
        err += {P.dT, dPdT_R};
    }

    // e(rho, P)
    {
        auto E = mix.energy_rP(z.rho, z.P, z.beta, {.deriv=true});

        double dEdR_P = derivative<1, 4>([&](double rho) -> double {
            return mix.energy_rP(rho, z.P, z.beta);
        }, z.rho, dR);
        double dEdP_R = derivative<1, 4>([&](double P) -> double {
            return mix.energy_rP(z.rho, P, z.beta);
        }, z.P, dP);

        err += {E.dR, dEdR_P};
        err += {E.dP, dEdP_R};
    }

    // e(rho, T)
    {
        auto E = mix.energy_rT(z.rho, z.T, z.beta, {.deriv=true});

        double dEdR_T = derivative<1, 4>([&](double rho) -> double {
            return mix.energy_rT(rho, z.T, z.beta);
        }, z.rho, dR);
        double dEdT_R = derivative<1, 4>([&](double T) -> double {
            return mix.energy_rT(z.rho, T, z.beta);
        }, z.T, dT);

        err += {E.dR, dEdR_T};
        err += {E.dT, dEdT_R};
    }

    // v(P, T)
    {
        auto V = mix.volume_PT(z.P, z.T, z.beta, {.deriv=true});

        double dVdP_T = derivative<1, 4>([&](double P) -> double {
            return mix.volume_PT(P, z.T, z.beta);
        }, z.P, dP);
        double dVdT_P = derivative<1, 4>([&](double T) -> double {
            return mix.volume_PT(z.P, T, z.beta);
        }, z.T, dT);

        err += {V.dP, dVdP_T};
        err += {V.dT, dVdT_P};
    }

    // e(P, T)
    {
        auto E = mix.energy_PT(z.P, z.T, z.beta, {.deriv=true});

        double dEdP_T = derivative<1, 4>([&](double P) -> double {
            return mix.energy_PT(P, z.T, z.beta);
        }, z.P, dP);
        double dEdT_P = derivative<1, 4>([&](double T) -> double {
            return mix.energy_PT(z.P, T, z.beta);
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
void test_sound_speed(const MixturePT& mix, State z) {
    dRdE P = mix.pressure_re(z.rho, z.e, z.beta, {.deriv=true});
    dRdP e = mix.energy_rP  (z.rho, z.P, z.beta, {.deriv=true});

    // Проверка скорости звука
    double c1 = std::sqrt(P.dR + P.val * P.dE / (z.rho * z.rho));
    double c2 = std::sqrt((P / (z.rho * z.rho) - e.dR) / e.dP);
    double c3 = mix.sound_speed_rP(z.rho, z.P, z.beta);
    double c4 = mix.sound_speed_re(z.rho, z.e, z.beta);

    ErrorList err = {
            {c1, c2},
            {c1, c3},
            {c1, c4}
    };

    if (err.is_ok(1.0e-12)) {
        std::cout << "\tSound Speed.  Max error: " << std::setw(12) << err.max() << ".\tOK!\n";
    } else {
        std::cout << "\tSound Speed.  Max error: " << std::setw(12) << err.max() << ".\tNOT OK!!\n";
        err.print("\t\t");
    }
}

// Аппроксимация двучленным УрС
void test_stiffened_gas(const MixturePT& mix, State z) {
    StiffenedGas sg = mix.stiffened_gas(z.rho, z.P, z.beta);

    dRdE P = mix.pressure_re(z.rho, z.e, z.beta, {.deriv=true});
    dRdP e = mix.energy_rP  (z.rho, z.P, z.beta, {.deriv=true});
    double c = mix.sound_speed_rP(z.rho, z.P, z.beta);

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

    if (err.is_ok(1.0e-12)) {
        std::cout << "\tStiffenedGas. Max error: " << std::setw(12) << err.max() << ".\tOK!\n";
    } else {
        std::cout << "\tStiffenedGas. Max error: " << std::setw(12) << err.max() << ".\tNOT OK!!\n";
        err.print("\t\t");
    }
}

/// @brief Полное тестирование
void test_mixture(MixturePT& mix, State z) {
    std::cout << "--------------------------------------------------------\n";
    z.print();
    std::cout << "--------------------------------------------------------\n";

    std::cout << std::scientific << std::setprecision(2);

    test_consistency(mix, z);
    test_derivatives(mix, z);
    test_sound_speed(mix, z);
    test_stiffened_gas(mix, z);

    std::cout << "--------------------------------------------------------\n";
}

// Смесь идеальных газов
MixturePT get_mixture_1(bool old_style) {
    MixturePT mixture(old_style);
    mixture += IdealGas::create(gam_1, Cv_1);
    mixture += IdealGas::create(gam_2, Cv_2);
    mixture += IdealGas::create(gam_3, Cv_3);

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "\tMaterial[0]: IdealGas(γ=" << gam_1 << ", cv=" << Cv_1 << ")\n";
    std::cout << "\tMaterial[1]: IdealGas(γ=" << gam_2 << ", cv=" << Cv_2 << ")\n";
    std::cout << "\tMaterial[2]: IdealGas(γ=" << gam_3 << ", cv=" << Cv_3 << ")\n";

    return mixture;
}

// Сложная смесь
MixturePT get_mixture_2(bool old_style) {
    MixturePT mixture(old_style);
    mixture += IdealGas::create("Air");
    mixture += StiffenedGas::create("Water2");
    mixture += StiffenedGas::create("Copper");

    std::cout << "\tMaterial[0]: IdealGas(\"Air\")\n";
    std::cout << "\tMaterial[1]: StiffenedGas(\"Water2\")\n";
    std::cout << "\tMaterial[2]: MieGruneisen(\"Cu\")\n";

    return mixture;
}

int main() {
    std::cout << "Mixture:\n";
    MixturePT mixture = get_mixture_2(false);
    //MixturePT mixture = get_mixture_2(true);

    // Параметры (не согласованы, и не обязаны)
    double rho_test{NAN}, e_test{NAN}, P_test{NAN}, T_test{NAN};
    Fractions alpha, beta;

    int mix_case = 0;
    switch (mix_case) {
        case 0: {
            // Простой тест, адекватная смесь
            rho_test = 2.13_g_cm3;
            e_test   = 0.29_MJ_kg;
            P_test   = 90.0_MPa;
            T_test   = 400.0_C;

            alpha = {0.5, 0.3, 0.2};
            beta  = {0.11, 0.07, 0.82};
            break;
        }
        case 1: {
            // Случай с давлением около нуля
            rho_test = 2.73e-5_g_cm3;
            e_test   = 2.71_MJ_kg;
            P_test   = 0.01_Pa;
            T_test   = 20.0_C;

            alpha = {0.9998, 0.0001, 0.3};
            beta  = {0.0004, 0.99, 0.0086};
            break;
        }
        case 2: {
            // Случай с нулевой объемной долей
            // при ненулевой массовой доле
            rho_test = 2.13_g_cm3;
            e_test   = 0.29_MJ_kg;
            P_test   = 90.0_MPa;
            T_test   = 400.0_C;

            alpha = {0.5, 0.3, 0.2};
            beta  = {0.11, 0.07, 0.82};
            break;
        }
        case 3: {
            // Точные значения для смеси идеальных газов
            P_test = 1.0_MPa;
            T_test = 300.0;

            beta = {0.1, 0.2, 0.7};
            double GC = beta[0] * (gam_1 - 1.0) * Cv_1 +
                        beta[1] * (gam_2 - 1.0) * Cv_2 +
                        beta[2] * (gam_3 - 1.0) * Cv_3;
            rho_test = P_test / (GC * T_test);
            e_test = (beta[0] * Cv_1 + beta[1] * Cv_2 + beta[2] * Cv_3) * T_test;

            ScalarSet rhos = {1.0 / ((gam_1 - 1.0) * Cv_1),
                              1.0 / ((gam_2 - 1.0) * Cv_2),
                              1.0 / ((gam_3 - 1.0) * Cv_3)};
            rhos *= P_test / T_test;

            alpha = {beta[0] * rho_test / rhos[0],
                     beta[1] * rho_test / rhos[1],
                     beta[2] * rho_test / rhos[2]};
            break;
        }
        default:
            throw std::runtime_error("Unknown test");
    }

    if (1) {
        std::cout << "\nTest mixture (explicit, (P, T))\n";
        State z = State::Explicit(mixture, alpha, P_test, T_test);
        test_mixture(mixture, z);
    }

    if (1) {
        std::cout << "\nTest mixture (implicit, (rho, P))\n";
        State z = State::Implicit_rP(mixture, beta, rho_test, P_test);
        test_mixture(mixture, z);
    }

    if (1) {
        std::cout << "\nTest mixture (implicit, (rho, T))\n";
        State z = State::Implicit_rT(mixture, beta, rho_test, T_test);
        test_mixture(mixture, z);
    }

    if (1) {
        std::cout << "\nTest mixture (implicit, (rho, e))\n";
        State z = State::Implicit_re(mixture, beta, rho_test, e_test);
        test_mixture(mixture, z);
    }

    return 0;
}