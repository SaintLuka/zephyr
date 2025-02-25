#include <iostream>

#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/matter/eos/stiffened_gas.h>
#include <zephyr/phys/matter/mixture_pt.h>

namespace zephyr::phys {

int MixturePT::size() const {
    return m_materials.size();
}

void MixturePT::clear() {
    m_materials.clear();
}

void MixturePT::append(Eos::Ref eos) {
    m_materials.emplace_back(eos);
}

void MixturePT::operator+=(Eos::Ref eos) {
    m_materials.emplace_back(eos);
}

// Обратный якобиан D(x, y) / D(T, P)
inline double inv_J(dPdT &x, dPdT &y) {
    return 1.0 / (x.dT * y.dP - x.dP * y.dT);
}

dRdE MixturePT::pressure_re(double rho, double e,
        const Fractions &beta, const MixOptions &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->pressure_re(rho, e, options[idx]);
    }

    return std::get<1>(find_rPT(rho, e, beta, options));
}

dRdT MixturePT::pressure_rT(double rho, double T,
        const Fractions &beta, const MixOptions &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->pressure_rT(rho, T, options[idx]);
    }

    return std::get<1>(find_rP_rT(rho, T, beta, options));
}

dRdT MixturePT::energy_rT(double rho, double T,
        const Fractions& beta, const MixOptions &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->energy_rT(rho, T, options[idx]);
    }

    return std::get<1>(find_reP_rT(rho, T, beta, options));
}

dRdP MixturePT::energy_rP(double rho, double P,
        const Fractions &beta, const MixOptions &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->energy_rP(rho, P);
    }

    return std::get<1>(find_reT_rP(rho, P, beta, options));
}

double MixturePT::temperature_rP(double rho, double P,
        const Fractions &beta, const MixOptions &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->temperature_rP(rho, P);
    }

    return std::get<1>(find_rT_rP(rho, P, beta, options));
}

double MixturePT::sound_speed_re(double rho, double e,
        const Fractions &beta, const MixOptions &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->sound_speed_re(rho, e);
    }

    MixOptions opts(options); opts.deriv = true;
    auto[rhos, P, T] = find_rPT(rho, e, beta, opts);
    double c2 = P.dR + P.val * P.dE / (rho * rho);
    return std::sqrt(c2);
}

double MixturePT::sound_speed_rP(double rho, double P,
        const Fractions &beta, const MixOptions &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->sound_speed_rP(rho, P);
    }

    auto opts = options; opts.deriv = true;
    auto[rhos, e, T] = find_reT_rP(rho, P, beta, opts);
    double c2 = (P / (rho * rho) - e.dR) / e.dP;
    return std::sqrt(c2);
}

dPdT MixturePT::volume_PT(double P, double T,
        const Fractions &beta, const MixOptions &options) const {

    dPdT v = {0.0, 0.0, 0.0};
    for (int i = 0; i < size(); ++i) {
        if (beta.has(i)) {
            // есть начальное приближение
            double rho_i0 = options.rhos ? (*options.rhos)[i] : NAN;

            auto v_i = m_materials[i]->volume_PT(P, T, {.deriv=true, .rho0=rho_i0});
            v.val += beta[i] * v_i.val;
            v.dP  += beta[i] * v_i.dP;
            v.dT  += beta[i] * v_i.dT;
        }
    }
    return v;
}

dPdT MixturePT::energy_PT(double P, double T,
        const Fractions &beta, const MixOptions &options) const {

    dPdT e = {0.0, 0.0, 0.0};
    for (int i = 0; i < size(); ++i) {
        if (beta.has(i)) {
            // есть начальное приближение
            double rho_i0 = options.rhos ? (*options.rhos)[i] : NAN;

            auto e_i = m_materials[i]->energy_PT(P, T, {.deriv=true, .rho0=rho_i0});
            e.val += beta[i] * e_i.val;
            e.dP  += beta[i] * e_i.dP;
            e.dT  += beta[i] * e_i.dT;
        }
    }
    return e;
}

StiffenedGas MixturePT::stiffened_gas(double rho, double P,
        const Fractions &beta, const MixOptions &options) const {
    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->stiffened_gas(rho, P);
    }

    auto opts = options; opts.deriv = true;

    auto[rhos, e, T] = find_reT_rP(rho, P, beta, opts);

    double gamma = 1.0 + 1.0 / (rho * e.dP);
    double e0 = e + rho * e.dR;
    double P0 = -(rho * e.dR + P * e.dP) / (gamma * e.dP);

    return StiffenedGas(gamma, P0, e0, NAN);
}

double MixturePT::min_pressure(const Fractions &beta) const {
    double P_min = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < size(); ++i) {
        if (beta.has(i)) {
            P_min = std::max(P_min, m_materials[i]->min_pressure());
        }
    }
    return P_min;
}

void MixturePT::adjust_cv(double rho_ref, double P_ref, double T_ref) {
    for (auto& mat: m_materials) {
        mat->adjust_cv(rho_ref, P_ref, T_ref);
    }
}

void MixturePT::adjust_T0(double rho_ref, double P_ref, double T_ref) {
    for (auto& mat: m_materials) {
        mat->adjust_T0(rho_ref, P_ref, T_ref);
    }
}

MixturePT::triplet_re MixturePT::get_rPT(double rho, double e,
        const Fractions& beta, const MixOptions& options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        ScalarSet rhos = ScalarSet::Pure(idx, rho);
        double P = m_materials[idx]->pressure_re(rho, e);
        double T = m_materials[idx]->temperature_rP(rho, P);
        return {rhos, P, T};
    }

    return find_rPT(rho, e, beta, options);
}

MixturePT::triplet_rP MixturePT::get_reT(double rho, double P,
        const Fractions& beta, const MixOptions& options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        ScalarSet rhos = ScalarSet::Pure(idx, rho);
        double e = m_materials[idx]->energy_rP(rho, P);
        double T = m_materials[idx]->temperature_rP(rho, P);
        return {rhos, e, T};
    }

    return find_reT_rP(rho, P, beta, options);
}

double MixturePT::sound_speed_PT(double P, double T,
        const Fractions &beta, const MixOptions &options) const {

    auto v = volume_PT(P, T, beta, {.deriv=true, .rhos=options.rhos});
    auto e = energy_PT(P, T, beta, {.deriv=true, .rhos=options.rhos});
    double inv_D = inv_J(v, e);

    double c2 = inv_D * std::pow(v, 2) * (e.dT + P * v.dT);
    return std::sqrt(c2);
}

ScalarSet MixturePT::init_densities(const Fractions &beta, const MixOptions &options) const {
    ScalarSet rhos = ScalarSet::NaN();
    if (options.rhos != nullptr) {
        for (int i = 0; i < size(); ++i) {
            if (beta.has(i)) {
                rhos[i] = (*options.rhos)[i];
            }
        }
    } else {
        for (int i = 0; i < size(); ++i) {
            if (beta.has(i)) {
                rhos[i] = m_materials[i]->density();
            }
        }
    }
    return rhos;
}

MixturePT::doublet_rT MixturePT::find_rP_rT(double density, double temperature,
        const Fractions& beta, const MixOptions& options) const {
    return m_old ? find_rP_rT_old(density, temperature, beta, options) :
                   find_rP_rT_new(density, temperature, beta, options);
}

MixturePT::triplet_rT MixturePT::find_reP_rT(double density, double temperature,
        const Fractions& beta, const MixOptions& options) const {
    return m_old ? find_reP_rT_old(density, temperature, beta, options) :
                   find_reP_rT_new(density, temperature, beta, options);
}

MixturePT::doublet_rP MixturePT::find_rT_rP(double density, double pressure,
        const Fractions& beta, const MixOptions& options) const {
    if (m_old) {
        return find_rT_rP_old(density, pressure, beta, options);
    } else {
        return find_rT_rP_new(density, pressure, beta, options);
    }
}

MixturePT::triplet_rP MixturePT::find_reT_rP(double density, double pressure,
        const Fractions& beta, const MixOptions& options) const {
    if (m_old) {
        return find_reT_rP_old(density, pressure, beta, options);
    } else {
        return find_reT_rP_new(density, pressure, beta, options);
    }
}

MixturePT::triplet_re MixturePT::find_rPT(double density, double energy,
        const Fractions& beta, const MixOptions& options) const {
    if (m_old) {
        return find_rPT_old(density, energy, beta, options);
    } else {
        return find_rPT_new(density, energy, beta, options);
    }
}

namespace {
const int max_counter = 30;
const double epsilon_v = 1.0e-12;
const double epsilon_T = 1.0e-12;
const double epsilon_P = 1.0e-12;
const double epsilon_e = 1.0e-12;
}

// Старая схема. Итерации по температуре. Не проверяет на чистоту
MixturePT::doublet_rT MixturePT::find_rP_rT_old(double rho, double T,
        const Fractions& beta, const MixOptions& options) const {
    // Решаем уравнение
    // v_mix = sum_i beta_i v_i(P, T)
    double v_mix = 1.0 / rho;
    double P = std::isnan(options.P0) ? 1.0e5 : options.P0;
    double P_min = min_pressure(beta);

    // Начальное приближение объемных долей
    ScalarSet rhos = init_densities(beta, options);

    //std::cout << "\tinitial \tP: " << P << ".\tϱ: " << rhos << "\n";

    dPdT v{NAN};
    int counter = 0;
    while (counter < max_counter) {
        double err = 0.0;
        v = {0.0, 0.0, 0.0};
        for (int i = 0; i < size(); ++i) {
            if (beta.has(i)) {
                auto v_i = m_materials[i]->volume_PT(P, T, {.deriv=true, .rho0=rhos[i]});
                v.val += beta[i] * v_i.val;
                v.dP  += beta[i] * v_i.dP;
                v.dT  += beta[i] * v_i.dT;

                err = std::max(err, std::abs(1.0 - rhos[i] * v_i));
                rhos[i] = 1.0 / v_i;
            }
        }
        err = std::max(err, std::abs(1.0 - rho * v));

        //std::cout << "\tcounter: " << counter << ".\tP: " << P << ".\tϱ: " << rhos << "; error: " << err << "\n";

        if (err < epsilon_v) {
            break;
        }

        // Итерация по Ньютону
        double dP = (v_mix - v) / v.dP;

        // Гарантирует выполнение P + dP > P_min
        if (P + dP < P_min) {
            // P + dP окажется на интервале (P_min, P)
            // равно (1 - a) * P_min + a * P

            const double a = 0.1;
            dP = (1.0 - a) * P_min + (a - 1.0) * P;
        }

        P += dP;

        ++counter;
    }

    dRdT P1 = P;
    if (options.deriv) {
        P1.dR = -std::pow(v.val, 2) / v.dP;
        P1.dT = -v.dT / v.dP;
    }

    return {rhos, P1};
}

MixturePT::triplet_rT MixturePT::find_reP_rT_old(double rho, double T,
        const Fractions& beta, const MixOptions& options) const {
    auto[rhos, P] = find_rP_rT_old(rho, T, beta, options);

    dPdT e1 = energy_PT(P, T, beta, {.deriv=options.deriv, .rho0=rho, .rhos=&rhos});
    dRdT e2 = e1.val;
    if (options.deriv) {
        e2.dR = e1.dP * P.dR;
        e2.dT = e1.dP * P.dT + e1.dT;
    }
    return {rhos, e2, P};
}

// Старая схема. Итерации по давлению. Не проверяет на чистоту
MixturePT::doublet_rP MixturePT::find_rT_rP_old(double rho, double P,
        const Fractions& beta, const MixOptions& options) const {
    // Решаем уравнение
    // v_mix = sum_i beta_i volume_i(P, T)
    double v_mix = 1.0 / rho;
    double T = std::isnan(options.T0) ? 300.0 : options.T0;

    // Начальное приближение плотностей
    ScalarSet rhos = init_densities(beta, options);

    //std::cout << "\tinitial: \tT: " << T << ".\tϱ: " << rhos << "\n";

    dPdT v{NAN};
    int counter = 0;
    while (counter < max_counter) {
        double err = 0.0;
        v = {0.0, 0.0, 0.0};
        for (int i = 0; i < size(); ++i) {
            if (beta.has(i)) {
                auto v_i = m_materials[i]->volume_PT(P, T, {.deriv=true, .rho0=rhos[i]});
                v.val += beta[i] * v_i.val;
                v.dP  += beta[i] * v_i.dP;
                v.dT  += beta[i] * v_i.dT;

                err = std::max(err, std::abs(1.0 - rhos[i] * v_i));
                rhos[i] = 1.0 / v_i;
            }
        }
        err = std::max(err, std::abs(1.0 - rho * v));

        //std::cout << "\tcounter: " << counter << ".\tT: " << T << ".\tϱ: " << rhos << "\n";

        if (err < epsilon_v) {
            break;
        }

        // Итерация по Ньютону
        T += (v_mix - v) / v.dT;

        ++counter;
    }

    dRdP T1 = T;
    if (options.deriv) {
        T1.dR = -std::pow(v.val, 2) / v.dT;
        T1.dP = -v.dP / v.dT;
    }

    return {rhos, T1};
}

MixturePT::triplet_rP MixturePT::find_reT_rP_old(double rho, double P,
        const Fractions& beta, const MixOptions& options) const {
    auto[rhos, T] = find_rT_rP_old(rho, P, beta, options);

    dPdT e1 = energy_PT(P, T, beta, {.deriv=options.deriv, .rho0=rho, .rhos=&rhos});
    dRdP e2 = e1.val;
    if (options.deriv) {
        e2.dR = e1.dT * T.dR;
        e2.dP = e1.dT * T.dP + e1.dP;
    }
    return {rhos, e2, T};
}

// Старая схема. Итерации по давлению и температуре. Не проверяет на чистоту
MixturePT::triplet_re MixturePT::find_rPT_old(double rho, double e_mix,
        const Fractions& beta, const MixOptions& options) const {
    // Решаем систему двух уравнений:
    //     v_mix = sum_i beta_i volume_i(P, T)
    //     e_mix = sum_i beta_i energy_i(P, T)

    double v_mix = 1.0 / rho;

    double P = std::isnan(options.P0) ? 1.0e5 : options.P0;
    double T = std::isnan(options.T0) ? 300.0 : options.T0;

    double P_min = min_pressure(beta);
    if (P < P_min) {
        P = P_min + 1;
    }

    // Начальное приближение объемных долей
    ScalarSet rhos = init_densities(beta, options);

    //std::cout << "\tinitial: \tP: " << P << ".\tT: " << T << ".\tϱ: " << rhos <<  "\n";

    int counter = 0;
    dPdT v{NAN}, e{NAN};
    while (counter < max_counter) {
        // Можно написать просто:
        //   v = volume_PT(P, T, beta, {.deriv=true, .rho0=rho, .rhos=&rhos});
        //   e = energy_PT(P, T, beta, {.deriv=true, .rho0=rho, .rhos=&rhos});
        // Но тогда будут двойные вычисления, и ещё rhos обновлять отдельно
        double err = 0.0;
        v = {0.0, 0.0, 0.0};
        e = {0.0, 0.0, 0.0};
        for (int i = 0; i < size(); ++i) {
            if (beta.has(i)) {
                auto v_i = m_materials[i]->volume_PT(P, T, {.deriv=true, .rho0=rhos[i]});
                v.val += beta[i] * v_i.val;
                v.dP  += beta[i] * v_i.dP;
                v.dT  += beta[i] * v_i.dT;

                err = std::max(err, std::abs(1.0 - rhos[i] * v_i));
                rhos[i] = 1.0 / v_i;

                auto e_i = m_materials[i]->energy_PT(P, T, {.deriv=true, .rho0=rhos[i]});
                e.val += beta[i] * e_i.val;
                e.dP  += beta[i] * e_i.dP;
                e.dT  += beta[i] * e_i.dT;
            }
        }
        err = std::max(err, std::abs(1.0 - rho * v));
        double err_e = std::abs(1.0 - e / e_mix);

        //std::cout << "\tcounter: " << counter << ".\tP: " << P << ".\tT: " << T << ".\tϱ: " << rhos <<  "\terror : " << err << "\n";

        if (err < epsilon_v && err_e < epsilon_e) {
            break;
        }

        double inv_D = inv_J(v, e);
        double dP = inv_D * ((e_mix - e) * v.dT - (v_mix - v) * e.dT);
        double dT = inv_D * ((v_mix - v) * e.dP - (e_mix - e) * v.dP);

        // Гарантирует выполнение P + dP > P_min
        double dP_lim = dP;
        if (P + dP < P_min) {
            const double a = 0.1;
            dP_lim = (1.0 - a) * P_min + (a - 1.0) * P;

            // P + dP_lim будет больше P_min
            // равно (1 - a) * P_min + a * P
        }
        double dT_lim = dP != 0.0 ? (dP_lim / dP) * dT : dT;

        P += dP_lim;
        T += dT_lim;

        ++counter;
    }

    dRdE P1 = P;
    if (options.deriv) {
        double inv_D = inv_J(v, e);
        P1.dE = inv_D * v.dT;
        P1.dR = inv_D * e.dT * std::pow(v, 2);
    }

    return {rhos, P1, T};
}


namespace Transform {
struct identity_map {
    inline double operator()(double x) const { return x; }
};
struct const_map {
    inline double operator()(double x) const { return 1.0; }
};
struct inverse_map {
    double val;
    inline double operator()(double x) const { return 1.0 / (x - val); }
};
struct inverse_inverse_map {
    double val;
    inline double operator()(double x) const { return 1.0 / x + val; }
};
struct inverse_deriv_map {
    double val;
    inline double operator()(double x) const { return -1.0 / std::pow(x + val, 2); }
};

std::tuple<identity_map, const_map, identity_map> Identity() {
    return { identity_map{}, const_map{}, identity_map{} };
}

std::tuple<inverse_map, inverse_deriv_map, inverse_inverse_map> Inverse(double P_min) {
    return { inverse_map{P_min}, inverse_deriv_map{P_min}, inverse_inverse_map{P_min} };
}
}

// Новая схема. Итерации по объемным долям. Не проверяет на чистоту
MixturePT::doublet_rT MixturePT::find_rP_rT_new(double rho, double T,
        const Fractions& beta, const MixOptions& options) const {
    // Решается система уравнений
    //   sum beta_i v_i = v,
    //   P_i (v_i, T) = P_j (v_j, T)

    // Число материалов в задаче
    int n_fractions = beta.count();

    double v = 1.0 / rho;

    double P_min = min_pressure(beta);
    //auto[f, df, inv_f] = Transform::Inverse(P_min);
    auto[f, df, inv_f] = Transform::Identity();

    // Начальное приближение плотностей (с нормировкой)
    ScalarSet rhos = init_densities(beta, options);

    //std::cout << "\tinitial: \tϱ: " << rhos << "\n";

    // Приращения объемных долей на итерациях
    ScalarSet dv;

    dRdT P = NAN;

    int counter = 0;
    while (counter < max_counter) {
        // Расчет вспомогательных величин
        double A_xi = 0.0;
        double B_xi = 0.0;
        double C_xi = 0.0;
        double v_hat = 0.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                double rho_i = rhos[i];
                dRdT P_i = m_materials[i]->pressure_rT(rho_i, T, {.deriv=true});

                double xi = beta[i] / (rho_i * rho_i * df(P_i) * P_i.dR);
                A_xi += xi;
                B_xi += xi * f(P_i);
                if (options.deriv) {
                    C_xi += xi * df(P_i) * P_i.dT;
                }

                v_hat += beta[i] / rho_i;
            }
        }

        // Вычисляем приращения dv
        // Находим константу chi: v + chi * dv > 0.
        double chi = 1.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                double rho_i = rhos[i];
                dRdT P_i = m_materials[i]->pressure_rT(rho_i, T, {.deriv=true});

                double coeff = 1.0 / (rho_i * rho_i * df(P_i) * P_i.dR);
                dv[i] = coeff * (f(P_i) - (B_xi + v_hat - v) / A_xi);

                if (dv[i] < 0.0) {
                    chi = std::min(chi, 0.95 / (rhos[i] * (-dv[i])));
                }
            }
        }

        // Поправляем приращения
        dv.arr() *= chi;

        // Вычисляем ошибку
        double err = 0.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                err = std::max(err, std::abs(rhos[i] * dv[i]));
                rhos[i] = rhos[i] / (1.0 + rhos[i] * dv[i]);
            }
        }
        err = std::max(err, std::abs(1.0 - rho * v_hat));

        //std::cout << "\tcounter: " << counter << ".\tP: " << inv_f(B_xi / A_xi) << ".\tϱ: " << rhos << "; \terror: " << err << "\n";

        if (err < epsilon_v) {
            P = inv_f(B_xi / A_xi);
            if (options.deriv) {
                P.dR = 1.0 / (rho * rho * df(P) * A_xi);
                P.dT = C_xi / (df(P) * A_xi);
            }
            break;
        }

        ++counter;
    }

    return {rhos, P};
}

MixturePT::triplet_rT MixturePT::find_reP_rT_new(double rho, double T,
        const Fractions& beta, const MixOptions& options) const {
    auto[rhos, P] = find_rP_rT_new(rho, T, beta, options);

    dPdT e1 = {0.0, 0.0, 0.0};
    for (int i = 0; i < size(); ++i) {
        if (beta.has(i)) {
            dRdT e_i = m_materials[i]->energy_rT(rhos[i], T, options[i]);
            e1.val += beta[i] * e_i.val;
            if (options.deriv) {
                dRdT P_i = m_materials[i]->pressure_rT(rhos[i], T, options[i]);

                e1.dP += beta[i] * e_i.dR / P_i.dR;
                e1.dT += beta[i] * (e_i.dT - e_i.dR * P_i.dT / P_i.dR);
            }
        }
    }

    dRdT e2 = e1.val;
    if (options.deriv) {
        e2.dR = e1.dP * P.dR;
        e2.dT = e1.dP * P.dT + e1.dT;
    }

    return {rhos, e2, P};
}

// Новая схема. Итерации по объемным долям и температуре. Не проверяет на чистоту
MixturePT::doublet_rP MixturePT::find_rT_rP_new(double rho, double P,
        const Fractions& beta, const MixOptions& options) const {
    // Решается система уравнений
    //   sum beta_i v_i = v,
    //   P_i (v_i, T) = P

    // Число материалов в задаче
    int n_fractions = beta.count();

    double v = 1.0 / rho;

    dRdP T = std::isnan(options.T0) ? 300.0 : options.T0;

    // Начальное приближение плотностей (с нормировкой)
    ScalarSet rhos = init_densities(beta, options);

    //std::cout << "\tinitial. \tT: " << T.val << ".\tϱ: " << rhos << "\n";

    // Приращения объемных долей на итерациях
    ScalarSet dv;

    int counter = 0;
    while (counter < max_counter) {
        // Расчет вспомогательных величин
        double A_xi = 0.0;
        double B_xi = 0.0;
        double C_xi = 0.0;
        double v_hat = 0.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                double rho_i = rhos[i];
                dRdT P_i = m_materials[i]->pressure_rT(rho_i, T, {.deriv=true});

                double xi = beta[i] / (rho_i * rho_i * P_i.dR);
                A_xi += xi;
                B_xi += xi * P_i;
                C_xi += xi * P_i.dT;

                v_hat += beta[i] / rho_i;
            }
        }

        double dT = (v - v_hat + A_xi * P - B_xi) / C_xi;

        // Вычисляем приращения dv
        // Находим константу chi: v + chi * dv > 0.
        double chi = 1.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                double rho_i = rhos[i];
                dRdT P_i = m_materials[i]->pressure_rT(rho_i, T, {.deriv=true});

                double coeff = 1.0 / (rho_i * rho_i * P_i.dR);
                dv[i] = coeff * (P_i - P + P_i.dT * dT);

                if (dv[i] < 0.0) {
                    chi = std::min(chi, 0.95 / (rhos[i] * (-dv[i])));
                }
            }
        }

        // 0 < T + chi * dT < 2T
        chi = std::min(chi, 0.95 * T / std::abs(dT));

        // Поправляем приращения
        dv *= chi;
        dT *= chi;

        // Вычисляем ошибку
        double err = 0.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                err = std::max(err, std::abs(rhos[i] * dv[i]));
                rhos[i] = rhos[i] / (1.0 + rhos[i] * dv[i]);
            }
        }
        double err_T = std::abs(dT / T);

        T.val += dT;

        //std::cout << "\tcounter: " << counter << ".\tT: " << T.val << ".\tϱ: " << rhos << "; error: " << err << "\n";

        if (err < epsilon_v && err_T < epsilon_T) {
            if (options.deriv) {
                T.dP = A_xi / C_xi;
                T.dR = -1.0 / (rho * rho * C_xi);
            }
            break;
        }

        ++counter;
    }

    return {rhos, T};
}

// Старая схема. Итерации по давлению. Не проверяет на чистоту
MixturePT::triplet_rP MixturePT::find_reT_rP_new(double rho, double P,
        const Fractions& beta, const MixOptions& options) const {
    auto[rhos, T] = find_rT_rP_new(rho, P, beta, options);

    dPdT e1 = {0.0, 0.0, 0.0};
    for (int i = 0; i < size(); ++i) {
        if (beta.has(i)) {
            dRdT e_i = m_materials[i]->energy_rT(rhos[i], T, options[i]);
            e1.val += beta[i] * e_i.val;
            if (options.deriv) {
                dRdT P_i = m_materials[i]->pressure_rT(rhos[i], T, options[i]);

                e1.dP += beta[i] * e_i.dR / P_i.dR;
                e1.dT += beta[i] * (e_i.dT - e_i.dR * P_i.dT / P_i.dR);
            }
        }
    }

    dRdP e2 = e1.val;
    if (options.deriv) {
        e2.dR = e1.dT * T.dR;
        e2.dP = e1.dT * T.dP + e1.dP;
    }
    return {rhos, e2, T};
}

// Новая схема. Итерации по объемным долям и температуре. Не проверяет на чистоту
MixturePT::triplet_re MixturePT::find_rPT_new(double rho, double e,
        const Fractions& beta, const MixOptions& options) const {
    // Решается система уравнений
    //   sum beta_i v_i = v,
    //   P_i (v_i, T) = P_j (v_j, T)
    //   sum beta_i e_i (v_i, T) = e

    dRdE P = NAN;
    dRdE T = std::isnan(options.T0) ? 300.0 : options.T0;

    double v = 1.0 / rho;

    // Число материалов в задаче
    int n_fractions = beta.count();

    // Начальное приближение объемных долей
    ScalarSet rhos = init_densities(beta, options);

    // Приращения объемных долей на итерациях
    ScalarSet dv;

    //std::cout << "\tinitial: \tT: " << T.val << ".\tϱ: " << rhos << "\n";

    int counter = 0;
    while (counter < max_counter) {
        // Расчет вспомогательных величин
        double A_xi{0.0}, A_eta{0.0};
        double B_xi{0.0}, B_eta{0.0};
        double C_xi{0.0}, C_eta{0.0};
        double v_hat = 0.0;
        double e_hat = 0.0;
        double e_hat_T = 0.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                double rho_i = rhos[i];
                dRdT P_i = m_materials[i]->pressure_rT(rho_i, T, {.deriv=true});
                dRdT e_i = m_materials[i]->energy_rT  (rho_i, T, {.deriv=true});

                double xi = beta[i] / (rho_i * rho_i * P_i.dR);
                double eta = beta[i] * e_i.dR / P_i.dR;

                A_xi += xi;
                B_xi += xi * P_i;
                C_xi += xi * P_i.dT;

                A_eta += eta;
                B_eta += eta * P_i;
                C_eta += eta * P_i.dT;

                v_hat += beta[i] / rho_i;
                e_hat += beta[i] * e_i.val;
                e_hat_T += beta[i] * e_i.dT;
            }
        }

        // Приращение dT не ограничиваем
        double dT = (A_xi * B_eta - A_eta * B_xi + A_xi * (e - e_hat) + A_eta * (v - v_hat)) /
                    (A_eta * C_xi - A_xi * C_eta + A_xi * e_hat_T);

        double chi = 1.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                double rho_i = rhos[i];
                dRdT P_i = m_materials[i]->pressure_rT(rho_i, T, {.deriv=true});

                double coeff = 1.0 / (rho_i * rho_i * P_i.dR);

                dv[i] = coeff * (P_i + (v - v_hat - B_xi) / A_xi + (P_i.dT - C_xi / A_xi) * dT);

                if (dv[i] < 0.0) {
                    chi = std::min(chi, 0.95 / (rhos[i] * (-dv[i])));
                }
            }
        }

        dv *= chi;
        dT *= chi;

        // Вычисляем ошибку
        double err = 0.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                err = std::max(err, std::abs(rhos[i] * dv[i]));
                rhos[i] = rhos[i] / (1.0 + rhos[i] * dv[i]);
            }
        }
        err = std::max(err, std::abs(1.0 - rho * v_hat));
        double err_e = std::abs(1.0 - e / e_hat);

        T.val += dT;

        //std::cout << "\tcounter: " << counter << ".\tP: " << B_xi / A_xi << ".\tT: " << T.val << ".\tϱ: " << rhos << "\terror : " << err << "\n";

        if (err < epsilon_v && err_e < epsilon_e) {
            P = B_xi / A_xi;
            if (options.deriv) {
                P.dR = (e_hat_T - C_eta) / (rho * rho * (C_xi * A_eta + A_xi * (e_hat_T - C_eta)));
                P.dE = C_xi / (C_xi * A_eta + A_xi * (e_hat_T - C_eta));
            }
            break;
        }

        ++counter;
    }

    return {rhos, P, T};
}

} // namespace zephyr::phys