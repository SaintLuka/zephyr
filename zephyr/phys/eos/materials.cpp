#include <iostream>
#include <zephyr/phys/eos/materials.h>

namespace zephyr::phys {

inline double sqr(double x) {
    return x * x;
}

int Materials::size() const {
    return m_materials.size();
}

void Materials::clear() {
    m_materials.clear();
}

void Materials::append(Eos::Ptr eos) {
    m_materials.emplace_back(eos);
}

void Materials::operator+=(Eos::Ptr eos) {
    m_materials.emplace_back(eos);
}

Eos &Materials::operator[](int idx) {
    return *m_materials[idx];
}

const Eos &Materials::operator[](int idx) const {
    return *m_materials[idx];
}

// Обратный якобиан D(x, y) / D(T, P)
inline double inv_J(dPdT &x, dPdT &y) {
    return 1.0 / (x.dT * y.dP - x.dP * y.dT);
}

dRdE Materials::pressure_re(double rho, double eps,
                            const Fractions &beta, const Options &options) const {
    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->pressure_re(rho, eps, options);
    }

    auto pair = find_PT(rho, eps, beta, options);

    dRdE P{pair.P};
    if (options.deriv) {
        auto v = volume_pt(pair.P, pair.T, beta);
        auto e = energy_pt(pair.P, pair.T, beta);
        double inv_D = inv_J(v, e);

        P.dE = inv_D * v.dT;
        P.dR = inv_D * e.dT * sqr(v);
    }
    return P;
}

double Materials::energy_rp(double rho, double P,
                            const Fractions &beta, const Options &options) const {
    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->energy_rp(rho, P, options);
    }

    double T = temperature_rp(rho, P, beta, options);
    return energy_pt(P, T, beta);
}


std::pair<double, double> Materials::temperature_energy_rp(double rho, double P, const Fractions &beta, const Options &options) const {
    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return {m_materials[idx]->temperature_rp(rho, P, options), m_materials[idx]->energy_rp(rho, P, options)};
    }

    double T = temperature_rp(rho, P, beta, options);
    return {T, energy_pt(P, T, beta)};
}

double Materials::sound_speed_re(double rho, double eps,
                                 const Fractions &beta, const Options &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->sound_speed_re(rho, eps, options);
    }

    auto pair = find_PT(rho, eps, beta, options);
    return sound_speed_pt(pair.P, pair.T, beta);
}

double Materials::sound_speed_rp(double rho, double P,
                                 const Fractions &beta, const Options &options) const {
    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->sound_speed_rp(rho, P, options);
    }

    double T = temperature_rp(rho, P, beta, options);
    return sound_speed_pt(P, T, beta);
}

double Materials::pressure_rt(double rho, double T,
                              const Fractions &beta, const Options &options) const {
    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->pressure_rt(rho, T, options);
    }

    // Решаем vol = sum_i beta_i vol_i(P, T)
    double vol = 1.0 / rho;
    double P = std::isnan(options.P0) ? 1.0e5 : options.P0;
    double P_min = min_pressure(beta);

    double err = 1.0;
    int counter = 0;
    while (err > 1.0e-10 && counter < 30) {
        auto v = volume_pt(P, T, beta);
        double dP = (vol - v) / v.dP;
        err = std::abs(dP / P);
        P += dP;
        if (P < P_min) {
            // TODO: Долго сходится при больших P0, оптимизировать
            // Среднее между P_min и предыдущим значением
            P = 0.5 * (P_min + P - dP);
        }
        ++counter;
    }
    return P;
}

double Materials::temperature_rp(double rho, double P,
                                 const Fractions &beta, const Options &options) const {
    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->temperature_rp(rho, P, options);
    }

    // Решаем vol = sum_i beta_i vol_i(P, T)
    double vol = 1.0 / rho;
    double T = std::isnan(options.T0) ? 300.0 : options.T0;

    double err = 1.0;
    int counter = 0;
    while (err > 1.0e-10 && counter < 30) {
        auto v = volume_pt(P, T, beta);
        double dT = (vol - v) / v.dT;
        err = std::abs(dT / T);
        T += dT;
        ++counter;
    }

    return T;
}

dPdT Materials::volume_pt(double P, double T, const Fractions &beta) const {
    dPdT vol = {0.0, 0.0, 0.0};
    for (int i = 0; i < size(); ++i) {
        if (beta.has(i)) {
            auto v_i = m_materials[i]->volume_pt(P, T, {.deriv = true});
            vol.val += beta[i] * v_i.val;
            vol.dP += beta[i] * v_i.dP;
            vol.dT += beta[i] * v_i.dT;
        }
    }
    return vol;
}

dPdT Materials::energy_pt(double P, double T, const Fractions &beta) const {
    dPdT eps = {0.0, 0.0, 0.0};
    for (int i = 0; i < size(); ++i) {
        if (beta.has(i)) {
            auto e_i = m_materials[i]->energy_pt(P, T, {.deriv = true});
            eps.val += beta[i] * e_i.val;
            eps.dP += beta[i] * e_i.dP;
            eps.dT += beta[i] * e_i.dT;
        }
    }
    return eps;
}

StiffenedGas Materials::stiffened_gas(double rho, double P,
                                      const Fractions &beta, const Options &options) const {
    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->stiffened_gas(rho, P, options);
    }

    double T = temperature_rp(rho, P, beta, options);
    auto vol = volume_pt(P, T, beta);
    auto eps = energy_pt(P, T, beta);
    double inv_D = inv_J(vol, eps);

    double gamma = 1.0 + vol * vol.dT * inv_D;
    double eps_0 = eps - vol * eps.dT / vol.dT;
    double P0 = (vol * eps.dT * inv_D - P) / gamma;

    return StiffenedGas(gamma, P0, eps_0, NAN);
}

double Materials::min_pressure(const Fractions &beta) const {
    double P_min = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < size(); ++i) {
        if (beta.has(i)) {
            P_min = std::max(P_min, m_materials[i]->min_pressure());
        }
    }
    return P_min;
}

PairPT Materials::find_PT(double rho, double eps,
                          const Fractions &beta, const Options &options) const {

    // Решаем vol = sum_i beta_i vol_i(P, T)
    //        eps = sum_i beta_i eps_i(P, T)

    double vol = 1.0 / rho;
    double P = std::isnan(options.P0) ? 1.0e5 : options.P0;
    double T = std::isnan(options.T0) ? temperature_rp(rho, P, beta) : options.T0;

    double P_min = min_pressure(beta);
    if (P < P_min)
        P = P_min + 1;

    double err = 1.0;
    int counter = 0;
    while (err > 1.0e-12 && counter < 30 && !std::isnan(P)) {
        auto v = volume_pt(P, T, beta);
        auto e = energy_pt(P, T, beta);

        double inv_D = inv_J(v, e);
        double dP = inv_D * ((eps - e) * v.dT - (vol - v) * e.dT);
        double dT = inv_D * ((vol - v) * e.dP - (eps - e) * v.dP);

        err = std::abs(dP / P) + std::abs(dT / T);

        if (P + dP < P_min) {
            // TODO: Долго сходится при больших P0, оптимизировать
            double dP_lim = 0.5 * (dP - (P - P_min) + sqrt(sqr(dP + P - P_min) + 0.01 * sqr(P - P_min)));
            double T_lim = dP_lim * dT / dP;

            P += dP_lim;
            T += T_lim;
        } else {
            P += dP;
            T += dT;
        }
        ++counter;
    }

    return {P, T};
}

double Materials::sound_speed_pt(double P, double T, const Fractions &beta) const {
    auto v = volume_pt(P, T, beta);
    auto e = energy_pt(P, T, beta);
    double inv_D = inv_J(v, e);

    double c2 = inv_D * sqr(v) * (e.dT + P * v.dT);
    return std::sqrt(c2);
}

}