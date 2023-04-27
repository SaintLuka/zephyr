#include <zephyr/math/solver/riemann.h>

namespace zephyr { namespace math {

using namespace zephyr::phys;
using namespace smf;

using cref = RiemannSolver::cref;

inline double sqr(double x) {
    return x * x;
}

/// @brief Скорость звука от плотности (rho) и давления (p)
inline double eos_sound(cref rho, cref p, cref gamma, cref p0) {
    return std::sqrt(gamma * (p + p0) / rho);
}

/// @brief Плотность от давления (p) и скорости звука (c)
inline double eos_density(cref p, cref c, cref gamma, cref p0) {
    return gamma * (p + p0) / sqr(c);
}

/// @brief Внутренняя энергия от плотности (rho) и давления (p)
inline double eos_energy(cref rho, cref p, cref gamma, cref p0, cref eps0) {
    return eps0 + (p + gamma * p0) / ((gamma - 1.0) * rho);
}

/// @brief Начальное давление для итераций по Ньютону
/// (акустическое приближение)
inline double p_init(
        cref rL, cref uL, cref pL, cref cL, cref p0L,
        cref rR, cref uR, cref pR, cref cR, cref p0R) {

    double rcL = rL * cL;
    double rcR = rR * cR;

    double P = (pL * rcR + pR * rcL + (uL - uR) * rcL * rcR) / (rcL + rcR);

    double Pmin = -std::min(p0L, p0R) + 1.0e-2 * std::min(pL + p0L, pR + p0R);

    return std::max(Pmin, P);
}

/// @brief Начальное давление для итераций по Ньютону
/// (две ударных волны)
inline double p_init_two_shocks(
        cref rL, cref uL, cref pL, cref cL, cref p0L,
        cref rR, cref uR, cref pR, cref cR, cref p0R) {

    return p_init(rL, uL, pL, cL, p0L,
                  rR, uR, pR, cR, p0R);
}

/// @brief Начальное давление для итераций по Ньютону
/// (две волны разрежения)
inline double p_init_two_rarefactions(
        cref rL, cref uL, cref pL, cref cL, cref p0L,
        cref rR, cref uR, cref pR, cref cR, cref p0R) {

    return p_init(rL, uL, pL, cL, p0L,
                  rR, uR, pR, cR, p0R);
}

inline double coeff_AK(cref rK, cref gK) {
    return 2.0 / ((gK + 1) * rK);
}

inline double coeff_BK(cref pK, cref gK, cref p0K) {
    return (gK - 1.0) * pK / (gK + 1.0) + 2.0 * gK * p0K / (gK + 1.0);
}

inline double coeff_GK(cref gamma) {
    return (gamma - 1.0) / (2.0 * gamma);
}

inline double func_fK(
        cref P, cref rK, cref pK, cref cK,
        cref gK, cref p0K, cref AK, cref BK, cref GK) {
    if (P >= pK) {
        return (P - pK) * std::sqrt(AK / (P + BK));
    } else {
        return 2.0 * cK / (gK - 1.0) * (std::pow((P + p0K) / (pK + p0K), GK) - 1.0);
    }
}

inline double deriv_fK(
        cref P, cref rK, cref pK, cref cK,
        cref gK, cref p0K, cref AK, cref BK, cref GK) {
    if (P >= pK) {
        return (1.0 - (P - pK) / (2.0 * (P + BK))) * std::sqrt(AK / (P + BK));
    } else {
        return std::pow((P + p0K) / (pK + p0K), GK - 1.0) / (rK * cK);
    }
}

/// @struct Минималистичная структура для хранения пары чисел:
/// давления P и скорости U на контактном разрыве при решении
/// задачи Римана о распаде разрыва
struct SolPU {
    double P;   ///< Давление на контакте
    double U;   ///< Скорость на контакте

    /// @brief Конструктор по умолчанию
    SolPU() : P(0.0 / 0.0), U(0.0 / 0.0) {}

    /// @brief Простейший конструктор
    SolPU(cref P, cref U) : P(P), U(U) {}
};

/// @brief Найти давление и скорость на контактном разрыве
/// @param rL, uL, pL Плотность, скорость, давление слева
/// @param rR, uR, pR Плотность, скорость, давление справа
/// @param cL, cR Скорость звука слева и справа
/// @param gL, p0L Параметра материала слева (двучленный УрС)
/// @param gR, p0R Параметры материала справа (двучленный УрС)
inline SolPU contact_p(
        cref rL, cref uL, cref pL, cref cL, cref gL, cref p0L,
        cref rR, cref uR, cref pR, cref cR, cref gR, cref p0R) {

    //std::cout << "Version 3\n";

    // Константы не меняются при итерациях, поэтому вычислим их сразу
    double AL = coeff_AK(rL, gL);
    double AR = coeff_AK(rR, gR);
    double BL = coeff_BK(pL, gL, p0L);
    double BR = coeff_BK(pR, gR, p0R);
    double GL = coeff_GK(gL);
    double GR = coeff_GK(gR);

    // Базовая классификация
    double p_dno = std::max(-p0L, -p0R);
    double p_min = std::min(pL, pR);
    double p_max = std::max(pL, pR);

    double f_dno = func_fK(p_dno, rL, pL, cL, gL, p0L, AL, BL, GL) +
                   func_fK(p_dno, rR, pR, cR, gR, p0R, AR, BR, GR) + uR - uL;

    // Вакуумный случай
    if (f_dno > 0.0) {
        std::cout << "Вакуум\n";
        return {p_dno, 0.0};
    }

    // Классифицируем случай и выбираем начальное приближение
    double P = 0.5*(pL + pR);

    double f_min = func_fK(p_min, rL, pL, cL, gL, p0L, AL, BL, GL) +
                   func_fK(p_min, rR, pR, cR, gR, p0R, AR, BR, GR) + uR - uL;

    // Две волны разрежения
    if (f_min > 0.0) {
        //std::cout << "Две волны разрежения\n";

        // Использовать приближение для двух волн разрежения
        P = p_init_two_rarefactions(rL, uL, pL, cL, p0L, rR, uR, pR, cR, p0R);

    } else {

        double f_max = func_fK(p_max, rL, pL, cL, gL, p0L, AL, BL, GL) +
                       func_fK(p_max, rR, pR, cR, gR, p0R, AR, BR, GR) + uR - uL;

        if (f_max < 0.0) {
            // Две ударных волны
            //std::cout << "Две ударных волны\n";

            // Использовать приближение двух ударных волн
            P = p_init_two_shocks(rL, uL, pL, cL, p0L, rR, uR, pR, cR, p0R);

        } else {
            //std::cout << "Ударная волна и волна разрежения\n";

            // Использовать акустическое приближение
            P = p_init(rL, uL, pL, cL, p0L, rR, uR, pR, cR, p0R);
        }
    }

    const int max_iterations = 50;

    double DP1 = 1.0e-5 * std::abs(pL - pR);
    double DP2 = 1.0e-12 * std::abs(std::max(pL + p0L, pR + p0R));
    double DP = std::max(DP1, DP2);

    //std::cout << "P init:" << P << "\n";

    double fL, fR, dfL, dfR;

    for (int counter = 0; counter < max_iterations; ++counter) {
        fL = func_fK(P, rL, pL, cL, gL, p0L, AL, BL, GL);
        fR = func_fK(P, rR, pR, cR, gR, p0R, AR, BR, GR);

        dfL = deriv_fK(P, rL, pL, cL, gL, p0L, AL, BL, GL);
        dfR = deriv_fK(P, rR, pR, cR, gR, p0R, AR, BR, GR);

        double Pn = P - (fL + fR + uR - uL) / (dfL + dfR);

        //std::cout << "Pn (" << counter << "): " << Pn << "\n";

        if (std::abs(Pn - P) < DP) {
            break;
        }

        P = Pn;
    }

    double U = 0.5 * (uL + uR + fR - fL);

    return {P, U};
}

/// @brief Скорость ударной волны
/// @param rK, uK, pK Плотность, скорость и давление перед фронтом
/// @param U, P Скорость и давление за фронтом ударной волны
inline double shock_wave_speed(cref rK, cref uK, cref pK, cref U, cref P) {
    return uK + (P - pK) / (rK * (U - uK));
}

/// @brief Плотность за фронтом ударной волны
/// @param rK, uK, pK Плотность, скорость и давление перед фронтом
/// @param U, P Скорость и давление за фронтом ударной волны
inline double shock_wave_density(cref rK, cref uK, cref pK, cref U, cref P) {
    return rK * (P - pK) / ((P - pK) - rK * sqr(U - uK));
}

/// @brief Скорость звука за левой волной разрежения
/// @param U Скорость за волной разрежения
/// @param uL, cL, gL Скорость, скорость звука и показатель адиабаты
/// в невозмущенной области
inline double rarefaction_sound_L(cref U, cref uL, cref cL, cref gL) {
    return cL + 0.5 * (gL - 1.0) * (uL - U);
}

/// @brief Скорость звука за правой волной разрежения
/// @param U Скорость за волной разрежения
/// @param uR, cR, gR Скорость, скорость звука и показатель адиабаты
/// в невозмущенной области
inline double rarefaction_sound_R(cref U, cref uR, cref cR, cref gR) {
    return cR - 0.5 * (gR - 1.0) * (uR - U);
}

/// @brief Скорость звука на левой волне разрежения
/// @param xi Автомодельная переменная xi = x/t
inline double sound_lfan(cref uL, cref cL, cref gL, cref xi = 0.0) {
    return ((gL - 1.0) * (uL - xi) + 2.0 * cL) / (gL + 1.0);
}

/// @brief Скорость звука на правой волне разрежения
/// @param xi Автомодельная переменная xi = x/t
inline double sound_rfan(cref uR, cref cR, cref gR, cref xi = 0.0) {
    return ((gR - 1.0) * (xi - uR) + 2.0 * cR) / (gR + 1.0);
}

/// @brief Плотность на левой волне разрежения
/// @param xi Автомодельная переменная xi = x/t
inline double density_lfan(cref rL, cref uL, cref cL, cref gL, cref xi = 0.0) {
    return rL * std::pow(2.0 / (gL + 1.0) + ((gL - 1.0) * (uL - xi)) / ((gL + 1.0) * cL), 2.0 / (gL - 1.0));
}

/// @brief Плотность на правой волне разрежения
/// @param xi Автомодельная переменная xi = x/t
inline double density_rfan(cref rR, cref uR, cref cR, cref gR, cref xi = 0.0) {
    return rR * std::pow(2.0 / (gR + 1.0) - ((gR - 1.0) * (uR - xi)) / ((gR + 1.0) * cR), 2.0 / (gR - 1.0));
}

/// @brief Скорость на левой волне разрежения
/// @param xi Автомодельная переменная xi = x/t
inline double velocity_lfan(cref uL, cref cL, cref gL, cref xi = 0.0) {
    return ((gL - 1.0) * uL + 2.0 * (cL + xi)) / (gL + 1.0);
}

/// @brief Скорость на правой волне разрежения
/// @param xi Автомодельная переменная xi = x/t
inline double velocity_rfan(cref uR, cref cR, cref gR, cref xi = 0.0) {
    return ((gR - 1.0) * uR + 2.0 * (xi - cR)) / (gR + 1.0);
}

/// @brief Давление на левой волне разрежения
/// @param xi Автомодельная переменная xi = x/t
inline double pressure_lfan(cref uL, cref pL, cref cL, cref gL, cref p0L, cref xi = 0.0) {
    return (pL + p0L) * std::pow(2.0 / (gL + 1.0) + ((gL - 1.0) * (uL - xi)) / ((gL + 1.0) * cL), 2.0 * gL / (gL - 1.0)) - p0L;
}

/// @brief Давление на правой волне разрежения
/// @param xi Автомодельная переменная xi = x/t
inline double pressure_rfan(cref uR, cref pR, cref cR, cref gR, cref p0R, cref xi = 0.0) {
    return (pR + p0R) * std::pow(2.0 / (gR + 1.0) - ((gR - 1.0) * (uR - xi)) / ((gR + 1.0) * cR), 2.0 * gR / (gR - 1.0)) - p0R;
}

RiemannSolver::Solution RiemannSolver::solve(
        const PState &zL, const PState &zR, const Eos&eos) {

    cref rL = zL.density;
    cref uL = zL.velocity.x();
    cref pL = zL.pressure;
    cref gL = eos.stiff_gamma(rL, pL);
    cref p0L = eos.stiff_p0(rL, pL);

    cref rR = zR.density;
    cref uR = zR.velocity.x();
    cref pR = zR.pressure;
    cref gR = eos.stiff_gamma(rR, pR);
    cref p0R = eos.stiff_p0(rR, pR);

    return RiemannSolver::solve(
            rL, uL, pL, gL, p0L,
            rR, uR, pR, gR, p0R);
}

RiemannSolver::Solution RiemannSolver::solve(
        cref rL, cref uL, cref pL, cref gL, cref p0L,
        cref rR, cref uR, cref pR, cref gR, cref p0R) {

    double cL = eos_sound(rL, pL, gL, p0L);
    double cR = eos_sound(rR, pR, gR, p0R);

    SolPU PU = contact_p(
            rL, uL, pL, cL, gL, p0L,
            rR, uR, pR, cR, gR, p0R);

    cref P = PU.P;
    cref U = PU.U;

    if (U > 0.0) {
        // Положительная скорость на контакте

        if (pL < P) {
            // Слева ударная волна (УВ)
            double D = shock_wave_speed(rL, uL, pL, U, P);

            if (D > 0.0) {
                // Положительная скорость УВ
                return {rL, uL, pL};
            } else {
                // Отрицательная скорость УВ
                return {shock_wave_density(rL, uL, pL, U, P), U, P};
            }
        } else {
            // Слева волна разрежения
            double DL1 = uL - cL;
            if (DL1 > 0.0) {
                // Волна разрежения полностью справа
                return {rL, uL, pL};
            } else {
                double cl = rarefaction_sound_L(U, uL, cL, gL);
                double DL2 = U - cl;
                if (DL2 < 0.0) {
                    // Волна разрежения полностью слева
                    return {eos_density(P, cl, gL, p0L), U, P};
                } else {
                    // Волна разрежения приходится на грань
                    return {
                            density_lfan(rL, uL, cL, gL),
                            velocity_lfan(uL, cL, gL),
                            pressure_lfan(uL, pL, cL, gL, p0L)
                    };
                }
            }
        }
    } else {
        // Отрицательная скорость на контакте

        if (P > pR) {
            // Справа ударная волна (УВ)
            double D = shock_wave_speed(rR, uR, pR, U, P);

            if (D < 0.0) {
                // Отрицательная скорость УВ
                return {rR, uR, pR};
            } else {
                // Положительная скорость УВ
                return {shock_wave_density(rR, uR, pR, U, P), U, P};
            }
        } else {
            // Справа волна разрежения

            double DR2 = uR + cR;

            if (DR2 < 0.0) {
                // Волна разрежения полностью слева
                return {rR, uR, pR};
            } else {
                double cr = rarefaction_sound_R(U, uR, cR, gR);
                double DR1 = U + cr;
                if (DR1 > 0.0) {
                    // Волна разрежения полностью справа
                    return {eos_density(P, cr, gR, p0R), U, P};
                } else {
                    // Волна разрежения приходится на грань
                    return {
                            density_rfan(rR, uR, cR, gR),
                            velocity_rfan(uR, cR, gR),
                            pressure_rfan(uR, pR, cR, gR, p0R)
                    };
                }
            }
        }
    }
}

RiemannSolver::RiemannSolver(
        const PState &zL, const PState &zR,
        const Eos &eos, double x_jump) :
        RiemannSolver(zL, zR, eos, eos, x_jump) { }

RiemannSolver::RiemannSolver(
        const PState &zL, const PState &zR,
        const Eos &eosL, const Eos &eosR,
        double x_jump) : x_jump(x_jump),
                         rL(zL.density), uL(zL.velocity.x()), pL(zL.pressure),
                         rR(zR.density), uR(zR.velocity.x()), pR(zR.pressure) {

    gL = eosL.stiff_gamma(rL, pL);
    p0L = eosL.stiff_p0(rL, pL);
    e0L = eosL.stiff_e0(rL, pL);

    gR = eosL.stiff_gamma(rR, pR);
    p0R = eosL.stiff_p0(rR, pR);
    e0R = eosL.stiff_e0(rR, pR);

    compute();
}

RiemannSolver::RiemannSolver(
        double rhoL, double uL, double pL, double gL, double p0L, double e0L,
        double rhoR, double uR, double pR, double gR, double p0R, double e0R,
        double x_jump) : x_jump(x_jump),
                         rL(rhoL), uL(uL), pL(pL), gL(gL), p0L(p0L), e0L(e0L),
                         rR(rhoR), uR(uR), pR(pR), gR(gR), p0R(p0R), e0R(e0R) {

    compute();
}

void RiemannSolver::compute() {
    cL = eos_sound(rL, pL, gL, p0L);
    cR = eos_sound(rR, pR, gR, p0R);

    auto UP = contact_p(
            rL, uL, pL, cL, gL, p0L,
            rR, uR, pR, cR, gR, p0R);

    P = UP.P;
    U = UP.U;

    // Скорости и плотности слева
    if (pL < P) {
        DL1 = DL2 = shock_wave_speed(rL, uL, pL, U, P);
        rl = shock_wave_density(rL, uL, pL, U, P);
    } else {
        double cl = rarefaction_sound_L(U, uL, cL, gL);
        DL1 = uL - cL;
        DL2 = U - cl;
        rl = eos_density(P, cl, gL, p0L);
    }

    // Скорости и плотности справа
    if (pR < P) {
        DR1 = DR2 = shock_wave_speed(rR, uR, pR, U, P);
        rr = shock_wave_density(rR, uR, pR, U, P);
    } else {
        double cr = rarefaction_sound_R(U, uR, cR, gR);
        DR1 = U + cr;
        DR2 = uR + cR;
        rr = eos_density(P, cr, gR, p0R);
    }
}

double RiemannSolver::sound_speed(double x, double t) const {
    if (t <= 0.0) {
        return x < x_jump ? cL : cR;
    }

    double xi = (x - x_jump) / t;
    if (xi < DL1) {
        return cL;
    }
    if (xi < DL2) {
        return sound_lfan(uL, cL, gL, xi);
    }
    if (xi < U) {
        return eos_sound(rl, P, gL, p0L);
    }
    if (xi < DR1) {
        return eos_sound(rr, P, gL, p0L);
    }
    if (xi < DR2) {
        return sound_rfan(uR, cR, gR, xi);
    }
    return cR;
}

double RiemannSolver::density(double x, double t) const {
    if (t <= 0.0) {
        return x < x_jump ? rL : rR;
    }

    double xi = (x - x_jump) / t;
    if (xi < DL1) {
        return rL;
    }
    if (xi < DL2) {
        return density_lfan(rL, uL, cL, gL, xi);
    }
    if (xi < U) {
        return rl;
    }
    if (xi < DR1) {
        return rr;
    }
    if (xi < DR2) {
        return density_rfan(rR, uR, cR, gR, xi);
    }
    return rR;
}

double RiemannSolver::velocity(double x, double t) const {
    if (t <= 0.0) {
        return x < x_jump ? uL : uR;
    }

    double xi = (x - x_jump) / t;
    if (xi < DL1) {
        return uL;
    }
    if (xi < DL2) {
        return velocity_lfan(uL, cL, gL, xi);
    }
    if (xi < DR1) {
        return U;
    }
    if (xi < DR2) {
        return velocity_rfan(uR, cR, gR, xi);
    }
    return uR;
}

double RiemannSolver::pressure(double x, double t) const {
    if (t <= 0.0) {
        return x < x_jump ? pL : pR;
    }

    double xi = (x - x_jump) / t;
    if (xi < DL1) {
        return pL;
    }
    if (xi < DL2) {
        return pressure_lfan(uL, pL, cL, gL, p0L, xi);
    }
    if (xi < DR1) {
        return P;
    }
    if (xi < DR2) {
        return pressure_rfan(uR, pR, cR, gR, p0R, xi);
    }
    return pR;
}

double RiemannSolver::energy(double x, double t) const {
    double g = x < x_jump + U * t ? gL : gR;
    double p0 = x < x_jump + U * t ? p0L : p0R;
    double e0 = x < x_jump + U * t ? e0L : e0R;

    double p = pressure(x, t);
    double rho = density(x, t);

    return eos_energy(rho, p, g, p0, e0);
}

} // namespace math
} // namespace zephyr