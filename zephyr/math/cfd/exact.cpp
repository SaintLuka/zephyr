#include <zephyr/math/cfd/exact.h>

namespace zephyr { namespace math {

using namespace zephyr::phys;
using namespace smf;

RiemannSolver::RiemannSolver(
        const PState &zL, const PState &zR,
        const Eos &eos) :
        RiemannSolver(zL, zR, eos, eos) {

    compute();
}

RiemannSolver::RiemannSolver(
        const PState &zL, const PState &zR,
        const Eos &eosL, const Eos &eosR) :
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
        double rhoR, double uR, double pR, double gR, double p0R, double e0R) :
        rL(rhoL), uL(uL), pL(pL), gL(gL), p0L(p0L), e0L(e0L),
        rR(rhoR), uR(uR), pR(pR), gR(gR), p0R(p0R), e0R(e0R) {

    compute();
}

using cref = RiemannSolver::cref;

inline double eos_sound(cref rho, cref p, cref gamma, cref p0) {
    return std::sqrt(gamma * (p + p0) / rho);
}

inline double const_gamma_p_1(cref gamma) {
    return 0.5 * (gamma + 1.0);
}

inline double const_gamma_m_1(cref gamma) {
    return 0.5 * (gamma - 1.0);
}

inline double const_kappa(cref gamma) {
    return (gamma - 1.0) / (2.0 * gamma);
}

inline double a_shock_g(cref P, cref rho, cref p, cref gamma, cref p0) {
    return std::sqrt(rho * (const_gamma_p_1(gamma) * (P + p0) + const_gamma_m_1(gamma) * (p + p0)));
}

inline double a_shock_gg(cref P, cref rho, cref p, cref g_p_1, cref g_m_1, cref p0) {
    return std::sqrt(rho * (g_p_1 * (P + p0) + g_m_1 * (p + p0)));
}

inline double a_rare_g(cref P, cref rho, cref p, cref c, cref gamma, cref p0) {
    double xi = (P + p0) / (p + p0);
    double kappa = const_kappa(gamma);
    return kappa * rho * c * (1.0 - xi) / (1.0 - std::pow(xi, kappa));
}

inline double a_rare_k(cref P, cref rho, cref p, cref c, cref kappa, cref p0) {
    double xi = (P + p0) / (p + p0);
    return kappa * rho * c * (1.0 - xi) / (1.0 - std::pow(xi, kappa));
}

inline double a_any(cref P, cref rho, cref p, cref c, cref gamma, cref p0) {
    if (P < p) {
        return a_rare_g(P, rho, p, c, gamma, p0);
    } else {
        return a_shock_g(P, rho, p, gamma, p0);
    }
}

inline double a_any(cref P, cref rho, cref p, cref c, cref kappa,
             cref gamma_p_1, cref gamma_m_1, cref p0) {
    if (P < p) {
        return a_rare_k(P, rho, p, c, kappa, p0);
    } else {
        return a_shock_gg(P, rho, p, gamma_p_1, gamma_m_1, p0);
    }
}

double RiemannSolver::contact_p_init(
        cref rL, cref uL, cref pL, cref cL,
        cref rR, cref uR, cref pR, cref cR) {

    double rcL = rL * cL;
    double rcR = rR * cR;

    return (pL * rcR + pR * rcL + (uL - uR) * rcL * rcR) / (rcL + rcR);
}

double RiemannSolver::contact_p_init(
        cref rL, cref uL, cref pL, cref gL, cref p0L,
        cref rR, cref uR, cref pR, cref gR, cref p0R) {

    double cL = eos_sound(rL, pL, gL, p0L);
    double cR = eos_sound(rR, pR, gR, p0R);

    return contact_p_init(rL, uL, pL, cL, rR, uR, pR, cR);
}

double RiemannSolver::contact_p(
        cref rL, cref uL, cref pL, cref cL, cref gL, cref p0L,
        cref rR, cref uR, cref pR, cref cR, cref gR, cref p0R) {

    const double eps = 1.0e-8;
    const int max_iterations = 30;

    double DP = eps * std::abs(pL - pR);

    double P = contact_p_init(rL, uL, pL, cL, rR, uR, pR, cR);

    double aL, aR, dP, Pk;

    double kL = const_kappa(gL);
    double gp1L = const_gamma_p_1(gL);
    double gm1L = const_gamma_m_1(gL);

    double kR = const_kappa(gR);
    double gp1R = const_gamma_p_1(gR);
    double gm1R = const_gamma_m_1(gR);

    for (int counter = 0; counter < max_iterations; ++counter) {
        aL = a_any(P, rL, pL, cL, kL, gp1L, gm1L, p0L);
        aR = a_any(P, rR, pR, cR, kR, gp1R, gm1R, p0R);

        Pk = (aR * pL + aL * pR + aL * aR * (uL - uR)) / (aL + aR);

        dP = std::abs(Pk - P);
        P = Pk;

        if (dP < DP) {
            break;
        }
    }

    return P;
}

double RiemannSolver::contact_p(
        cref rL, cref uL, cref pL, cref gL, cref p0L,
        cref rR, cref uR, cref pR, cref gR, cref p0R) {

    double cL = eos_sound(rL, pL, gL, p0L);
    double cR = eos_sound(rR, pR, gR, p0R);

    return contact_p(rL, uL, pL, cL, gL, p0L,
                     rR, uR, pR, cR, gR, p0R);
}

double RiemannSolver::contact_u(cref uL, cref pL, cref aL, cref uR, cref pR, cref aR) {
    return (aL * uL + aR * uR + pL - pR) / (aL + aR);
}

void RiemannSolver::compute() {
    cL = eos_sound(rL, pL, gL, p0L);
    cR = eos_sound(rR, pR, gR, p0R);

    P = contact_p(
            rL, uL, pL, cL, gL, p0L,
            rR, uR, pR, cR, gR, p0R);

    double aL = a_any(P, rL, pL, cL, gL, p0L);
    double aR = a_any(P, rR, pR, cR, gR, p0R);

    U = contact_u(uL, pL, aL, uR, pR, aR);

    // Скорсти и плотности слева
    if (pL < P) {
        DL1 = DL2 = uL - aL / rL;
        rl = rL * aL / (aL - rL * (uL - U));
    } else {
        DL1 = uL - cL;
        cl = cL + 0.5 * (gL - 1) * (uL - U);
        DL2 = U - cl;
        rl = gL * (P + p0L) / (cl * cl);
    }

    // Скорости и плотности справа
    if (pR < P) {
        DR1 = DR2 = uR + aR / rR;
        rr = rR * aR / (aR + rR * (uR - U));
    } else {
        cr = cR - 0.5 * (gR - 1) * (uR - U);
        DR1 = U + cr;
        rr = gL * (P + p0R) / (cr * cr);
        DR2 = uR + cR;
    }
}

double RiemannSolver::sound_speed(double x, double t) const {
    if (t <= 0.0) {
        return x < 0.0 ? cL : cR;
    }

    double xi = x / t;
    if (xi < DL1) {
        return cL;
    }
    if (xi < DL2) {
        return ((gL - 1.0) * (uL - xi) + 2.0 * cL) / (gL + 1.0);
    }
    if (xi < U) {
        return cl;
    }
    if (xi < DR1) {
        return cr;
    }
    if (xi < DR2) {
        return ((gR - 1.0) * (xi - uR) + 2.0 * cR) / (gR + 1.0);
    }
    return cR;
}

double RiemannSolver::density(double x, double t) const {
    if (t <= 0.0) {
        return x < 0.0 ? rL : rR;
    }

    double xi = x / t;
    if (xi < DL1) {
        return rL;
    }
    if (xi < DL2) {
        double c = sound_speed(x, t);
        double c2 = c * c;
        return std::pow(c2 * std::pow(rL, gL) / (gL * (pL + p0L)), 1.0 / (gL - 1.0));
    }
    if (xi < U) {
        return rl;
    }
    if (xi < DR1) {
        return rr;
    }
    if (xi < DR2) {
        double c = sound_speed(x, t);
        double c2 = c * c;
        return std::pow(c2 * std::pow(rR, gR) / (gR * (pR + p0R)), 1.0 / (gR - 1.0));
    }
    return rR;
}

double RiemannSolver::velocity(double x, double t) const {
    if (t <= 0.0) {
        return x < 0.0 ? uL : uR;
    }

    double xi = x / t;
    if (xi < DL1) {
        return uL;
    }
    if (xi < DL2) {
        return ((gL - 1.0) * uL + 2.0 * (xi + cL)) / (gL + 1.0);
    }
    if (xi < DR1) {
        return U;
    }
    if (xi < DR2) {
        return ((gR - 1.0) * uR + 2.0 * (xi - cR)) / (gR + 1.0);
    }
    return uR;
}

double RiemannSolver::pressure(double x, double t) const {
    if (t <= 0.0) {
        return x < 0.0 ? pL : pR;
    }

    double xi = x / t;
    if (xi < DL1) {
        return pL;
    }
    if (xi < DL2) {
        double rho = density(x, t);
        return (pL + p0L) * std::pow(rho / rL, gL) + p0L;
    }
    if (xi < DR1) {
        return P;
    }
    if (xi < DR2) {
        double rho = density(x, t);
        return (pR + p0R) * std::pow(rho / rR, gR) + p0R;
    }
    return pR;
}

double RiemannSolver::energy(double x, double t) const {
    double g = x < U * t ? gL : gR;
    double p0 = x < U * t ? p0L : p0R;
    double e0 = x < U * t ? e0L : e0R;

    double p = pressure(x, t);
    double rho = density(x, t);

    return e0 + (p + g * p0) / ((g - 1.0) * rho);
}

} // namespace math
} // namespace zephyr