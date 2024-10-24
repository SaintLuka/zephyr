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

dRdE Materials::pressure_re(double rho, double eps, const Fractions &beta,
        const Options &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        if (options.alpha) {
            options.alpha->set_pure(idx);
        }
        return m_materials[idx]->pressure_re(rho, eps, options);
    }

    auto pair = find_PT(rho, eps, beta, options);

    dRdE P{pair.P};
    if (options.deriv) {
        auto v = volume_PT(pair.P, pair.T, beta, {.deriv=true, .rho0=rho, .alpha=options.alpha});
        auto e = energy_PT(pair.P, pair.T, beta, {.deriv=true, .rho0=rho, .alpha=options.alpha});
        double inv_D = inv_J(v, e);

        P.dE = inv_D * v.dT;
        P.dR = inv_D * e.dT * sqr(v);
    }
    return P;
}

dRdE Materials::pressure_re2(double rho, double eps, const Fractions &beta,
                            const Options &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        if (options.alpha) {
            options.alpha->set_pure(idx);
        }
        return m_materials[idx]->pressure_re(rho, eps, options);
    }

    auto pair = find_PT_ver2(rho, eps, beta, options);

    dRdE P{pair.P};
    if (options.deriv) {
        auto v = volume_PT(pair.P, pair.T, beta, {.deriv=true, .rho0=rho, .alpha=options.alpha});
        auto e = energy_PT(pair.P, pair.T, beta, {.deriv=true, .rho0=rho, .alpha=options.alpha});
        double inv_D = inv_J(v, e);

        P.dE = inv_D * v.dT;
        P.dR = inv_D * e.dT * sqr(v);
    }
    return P;
}

double Materials::energy_rP(double rho, double P, const Fractions &beta,
        const Options &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        if (options.alpha) {
            options.alpha->set_pure(idx);
        }
        return m_materials[idx]->energy_rP(rho, P, options);
    }

    double T = temperature_rP(rho, P, beta, options);
    return energy_PT(P, T, beta, {.rho0=rho, .alpha=options.alpha});
}

double Materials::sound_speed_re(double rho, double eps, const Fractions &beta,
        const Options &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        if (options.alpha) {
            options.alpha->set_pure(idx);
        }
        return m_materials[idx]->sound_speed_re(rho, eps, options);
    }

    auto pair = find_PT(rho, eps, beta, options);
    return sound_speed_PT(pair.P, pair.T, beta, {.rho0=rho, .alpha=options.alpha});
}

double Materials::sound_speed_rP(double rho, double P, const Fractions &beta,
        const Options &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        if (options.alpha) {
            options.alpha->set_pure(idx);
        }
        return m_materials[idx]->sound_speed_rP(rho, P, options);
    }

    double T = temperature_rP(rho, P, beta, options);
    return sound_speed_PT(P, T, beta, {.rho0=rho, .alpha=options.alpha});
}

double Materials::pressure_rT_ver1(double rho, double T, const Fractions &beta,
        const Options &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        if (options.alpha) {
            options.alpha->set_pure(idx);
        }
        return m_materials[idx]->pressure_rT(rho, T, options);
    }

    // Решаем vol = sum_i beta_i vol_i(P, T)
    double vol = 1.0 / rho;
    double P = std::isnan(options.P0) ? 1.0e5 : options.P0;
    double P_min = min_pressure(beta);

    // Начальное приближение объемных долей
    Fractions alpha = options.alpha ? *options.alpha : beta;

    double err = 1.0;
    int counter = 0;
    while (err > 1.0e-12 && counter < 30) {
        auto v = volume_PT(P, T, beta, {.deriv=true, .rho0=rho, .alpha=&alpha});

        double dP = (vol - v) / v.dP;

        // Гарантирует выполнение P + dP > P_min
        if (P + dP < P_min) {
            // P + dP окажется на интервале (P_min, P)
            // равно (1 - a) * P_min + a * P

            const double a = 0.1;
            dP = (1.0 - a) * P_min + (a - 1.0) * P;
        }

        err = std::abs(dP / P);

        P += dP;

        // Обновить объемные доли
        update_alpha(rho, P, T, beta, alpha);

        //std::cout << "\titer " << counter << ".\tP: " << P << ".\ta: " << alpha << "\n";

        ++counter;
    }

    // Перенести объемные доли
    if (options.alpha) {
        *options.alpha = alpha;
    }

    return P;
}

double Materials::pressure_rT_ver2(double rho, double T, const Fractions &beta,
                              const Options &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        if (options.alpha) {
            options.alpha->set_pure(idx);
        }
        return m_materials[idx]->pressure_rT(rho, T, options);
    }

    // Решается система уравнений
    //   sum alpha_i = 1,
    //   P_i (rho_i, T) = P_j (rho_j, T),
    // где rho_i = beta_i rho / alpha_i

    // Число материалов в задаче
    int n_fractions = beta.count();

    // Начальное приближение объемных долей
    Fractions alpha = options.alpha ? *options.alpha : beta;

    // Приращения объемных долей на итерациях
    ScalarSet delta;

    double err = 1.0;
    int counter = 0;
    while (err > 1.0e-12 && counter < 30) {
        // Расчет вспомогательных величин
        double A_xi = 0.0;
        double B_xi = 0.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                double rho_i = (beta[i] / alpha[i]) * rho;
                dRdT P_i = m_materials[i]->pressure_rT(rho_i, T, {.deriv=true});

                double xi = alpha[i] / (rho_i * P_i.dR);
                A_xi += xi;
                B_xi += xi * P_i;
            }
        }

        // Вычисляем приращения alpha.
        // Сумма delta равна нулю по построению, а значит
        // сумма alpha будет равна единице.
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                double rho_i = (beta[i] / alpha[i]) * rho;
                dRdT P_i = m_materials[i]->pressure_rT(rho_i, T, {.deriv=true});

                double xi = alpha[i] / (rho_i * P_i.dR);
                delta[i] = xi * (P_i - B_xi / A_xi);
            }
        }

        // Поправляем условия, чтобы 0 < alpha_i + da_i < 1.
        // Найдем константу chi, потом умножим на неё.
        double chi = 1.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                if (alpha[i] + delta[i] > 1.0) {
                    chi = std::min(chi, 0.95 * (1.0 - alpha[i]) / (alpha[i] + delta[i] - 1.0));
                }
                else if (alpha[i] + delta[i] < 0.0) {
                    chi = std::min(chi, 0.95 * (alpha[i]) / (-delta[i]));
                }
            }
        }
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                delta[i] *= chi;
            }
        }

        err = 0.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                alpha[i] += delta[i];
                err = std::max(err, std::abs(delta[i]));

                assert(alpha[i] < 1.0);
                assert(alpha[i] > 0.0);
            }
        }

        ++counter;
    }

    // Просто усредняем. На самом деле, если метод сошелся,
    // тогда все значения P_i равны
    double P = 0.0;
    for (int i = 0; i < n_fractions; ++i) {
        if (beta.has(i)) {
            double rho_i = (beta[i] / alpha[i]) * rho;
            P += m_materials[i]->pressure_rT(rho_i, T);
        }
    }
    P /= beta.nonzero();

    // Перенести объемные доли
    if (options.alpha) {
        *options.alpha = alpha;
    }

    return P;
}

double Materials::pressure_rT(double rho, double T, const Fractions &beta,
                              const Options &options) const {

    return pressure_rT_ver2(rho, T, beta, options);
}

double Materials::temperature_rP(double rho, double P, const Fractions &beta,
        const Options &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        if (options.alpha) {
            options.alpha->set_pure(idx);
        }
        return m_materials[idx]->temperature_rP(rho, P, options);
    }

    // Решаем vol = sum_i beta_i vol_i(P, T)
    double vol = 1.0 / rho;
    double T = std::isnan(options.T0) ? 300.0 : options.T0;

    // Начальное приближение объемных долей
    Fractions alpha = options.alpha ? *options.alpha : beta;

    double err = 1.0;
    int counter = 0;
    while (err > 1.0e-10 && counter < 30) {
        auto v = volume_PT(P, T, beta, {.deriv=true, .rho0=rho, .alpha=&alpha});

        double dT = (vol - v) / v.dT;

        err = std::abs(dT / T);

        T += dT;

        // Обновить объемные доли
        update_alpha(rho, P, T, beta, alpha);

        ++counter;
    }

    // Перенести объемные доли
    if (options.alpha) {
        *options.alpha = alpha;
    }

    return T;
}

dPdT Materials::volume_PT(double P, double T, const Fractions &beta,
        const Options &options) const {

    dPdT vol = {0.0, 0.0, 0.0};
    for (int i = 0; i < size(); ++i) {
        if (beta.has(i)) {
            // есть начальное приближение
            double rho_i0 = NAN;
            if (!std::isnan(options.rho0) && options.alpha) {
                rho_i0 = beta[i] * options.rho0 / (*options.alpha)[i];
            }

            auto v_i = m_materials[i]->volume_PT(P, T, {.deriv = true, .rho0 = rho_i0});
            vol.val += beta[i] * v_i.val;
            vol.dP  += beta[i] * v_i.dP;
            vol.dT  += beta[i] * v_i.dT;
        }
    }
    return vol;
}

dPdT Materials::energy_PT(double P, double T, const Fractions &beta,
        const Options &options) const {

    dPdT eps = {0.0, 0.0, 0.0};
    for (int i = 0; i < size(); ++i) {
        if (beta.has(i)) {
            // есть начальное приближение
            double rho_i0 = NAN;
            if (!std::isnan(options.rho0) && options.alpha) {
                rho_i0 = beta[i] * options.rho0 / (*options.alpha)[i];
            }

            auto e_i = m_materials[i]->energy_PT(P, T, {.deriv = true, .rho0 = rho_i0});
            eps.val += beta[i] * e_i.val;
            eps.dP  += beta[i] * e_i.dP;
            eps.dT  += beta[i] * e_i.dT;
        }
    }
    return eps;
}

StiffenedGas Materials::stiffened_gas(double rho, double P,
        const Fractions &beta, const Options &options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        if (options.alpha) {
            options.alpha->set_pure(idx);
        }
        return m_materials[idx]->stiffened_gas(rho, P, options);
    }

    double T = temperature_rP(rho, P, beta, options);
    auto vol = volume_PT(P, T, beta, {.deriv=true, .rho0=rho, .alpha=options.alpha});
    auto eps = energy_PT(P, T, beta, {.deriv=true, .rho0=rho, .alpha=options.alpha});
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

void Materials::update_alpha(double rho, double P, double T,
                             const Fractions& beta, Fractions& alpha) const {

    // Обновить объемные доли
    for (int i = 0; i < size(); ++i) {
        if (beta[i] > 0.0) {
            // Можно оценить удельный объем компоненты
            double rho_i0 = NAN;
            if (alpha[i] > 0.0) {
                rho_i0 = beta[i] * rho / alpha[i];
            }
            double v_i0 = m_materials[i]->volume_PT(P, T, {.rho0 = rho_i0});

            alpha[i] = beta[i] * rho * v_i0;
        }
        else {
            alpha[i] = 0.0;
        }
    }
    alpha.normalize();
}

PairET Materials::find_eT(double rho, double P, const Fractions& beta,
        const Options& options) const {

    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        if (options.alpha) {
            options.alpha->set_pure(idx);
        }
        return {.e = m_materials[idx]->energy_rP(rho, P, options),
                .T = m_materials[idx]->temperature_rP(rho, P, options)};
    }

    double T = temperature_rP(rho, P, beta, {.deriv=false, .T0=options.T0, .alpha=options.alpha});
    double e = energy_PT(P, T, beta, {.deriv=false, .rho0=rho, .alpha=options.alpha});
    return {.e = e, .T = T};
}

PairPT Materials::find_PT_ver1(double rho, double eps, const Fractions &beta,
        const Options &options) const {

    // Решаем vol = sum_i beta_i vol_i(P, T)
    //        eps = sum_i beta_i eps_i(P, T)

    double vol = 1.0 / rho;
    double P = std::isnan(options.P0) ? 1.0e5 : options.P0;
    double T = std::isnan(options.T0) ? 300.0 : options.T0;

    double P_min = min_pressure(beta);
    if (P < P_min) {
        P = P_min + 1;
    }

    // Начальное приближение объемных долей
    Fractions alpha = options.alpha ? *options.alpha : beta;

    double err = 1.0;
    int counter = 0;
    while (err > 1.0e-12 && counter < 30 && !std::isnan(P)) {
        auto v = volume_PT(P, T, beta, {.deriv=true, .rho0=rho, .alpha=&alpha});
        auto e = energy_PT(P, T, beta, {.deriv=true, .rho0=rho, .alpha=&alpha});

        double inv_D = inv_J(v, e);
        double dP = inv_D * ((eps - e) * v.dT - (vol - v) * e.dT);
        double dT = inv_D * ((vol - v) * e.dP - (eps - e) * v.dP);

        err = std::abs(dP / P) + std::abs(dT / T);

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

        // Обновить объемные доли
        update_alpha(rho, P, T, beta, alpha);

        //std::cout << "iter: " << counter << "; P: " << P << "; T: " << T << "; a: " << alpha << "\n";

        ++counter;
    }

    // Перенести объемные доли
    if (options.alpha) {
        *options.alpha = alpha;
    }

    return {P, T};
}

PairPT Materials::find_PT_ver2(double rho, double e, const Fractions &beta,
                               const Options &options) const {
    using Eigen::VectorXd;
    using Eigen::MatrixXd;

    // Случай одного материала
    int cln = beta.index();
    if (cln >= 0) {
        if (options.alpha) {
            options.alpha->set_pure(cln);
        }
        double P = m_materials[cln]->pressure_re(rho, e, options);
        double T = m_materials[cln]->temperature_rP(rho, P, options);
        return {.P = P, .T = T};
    }

    double T = std::isnan(options.T0) ? 300.0 : options.T0;

    // Число материалов в задаче
    int n_fractions = beta.count();

    // Начальное приближение объемных долей
    Fractions alpha = options.alpha ? *options.alpha : beta;
    for (int i = 0; i < n_fractions; ++i) {
        if (beta.has(i) && alpha[i] == 0.0) {
            alpha[i] = 0.0001; // добавить материал
        }
    }
    alpha.normalize();

    // Приращения объемных долей на итерациях
    ScalarSet delta;

    // Индексы ненулевых материалов
    std::vector<int> idx;
    for (int i = 0; i < n_fractions; ++i) {
        if (beta.has(i)) {
            idx.push_back(i);
        }
    }

    // Количество ненулевых
    int nonzero = int(idx.size());

    double err = 1.0;
    int counter = 0;
    while (err > 1.0e-12 && counter < 30) {
        // Заполняем матрицу и правые части
        // TODO: Переписать с явными формулами?
        VectorXd F = VectorXd::Zero(nonzero + 1);
        MatrixXd A = MatrixXd::Zero(nonzero + 1, nonzero + 1);

        for (int m1 = 0; m1 < nonzero - 1; ++m1) {
            int m2 = m1 + 1;

            int i = idx[m1];
            int j = idx[m2];

            double rho_i = beta[i] * rho / alpha[i];
            double rho_j = beta[j] * rho / alpha[j];

            dRdT P_i = m_materials[i]->pressure_rT(rho_i, T, {.deriv=true});
            dRdT P_j = m_materials[j]->pressure_rT(rho_j, T, {.deriv=true});

            A(m1, m1) = + rho_i * P_i.dR / alpha[i];
            A(m1, m2) = - rho_j * P_j.dR / alpha[j];
            A(m1, nonzero) = P_j.dT - P_i.dT;

            F(m1) = P_i.val - P_j.val;
        }

        for (int m1 = 0; m1 < nonzero; ++m1) {
            A(nonzero - 1, m1) = 1.0;
        }
        F(nonzero - 1) = 0.0;

        F(nonzero) = e;
        for (int m1 = 0; m1 < nonzero; ++m1) {

            int i = idx[m1];

            double rho_i = beta[i] * rho / alpha[i];
            dRdT e_i = m_materials[i]->energy_rT(rho_i, T, {.deriv=true});

            A(nonzero, m1) -= beta[i] * rho_i * e_i.dR / alpha[i];
            A(nonzero, nonzero) += beta[i] * e_i.dT;

            F(nonzero) -= beta[i] * e_i.val;
        }

        // Решение системки на приращения
        VectorXd x = A.fullPivLu().solve(F);

        // Переносим приращения
        for (int i = 0; i < nonzero; ++i) {
            delta[idx[i]] = x[i];
        }
        double dT = x[nonzero];

        // Поправляем условия, чтобы 0 < alpha_i + da_i < 1.
        // Найдем константу chi, потом умножим на неё.
        double chi = 1.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                if (alpha[i] + delta[i] > 1.0) {
                    chi = std::min(chi, 0.95 * (1.0 - alpha[i]) / (delta[i]));
                }
                else if (alpha[i] + delta[i] < 0.0) {
                    chi = std::min(chi, 0.95 * (alpha[i]) / (-delta[i]));
                }

            }
        }
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                delta[i] *= chi;
            }
        }
        dT *= chi;

        err = 0.0;
        for (int i = 0; i < n_fractions; ++i) {
            if (beta.has(i)) {
                alpha[i] += delta[i];
                err = std::max(err, std::abs(delta[i]));

                assert(alpha[i] < 1.0);
                assert(alpha[i] > 0.0);
            }
        }
        T += dT;

        ++counter;
    }

    // Просто усредняем. На самом деле, если метод сошелся,
    // тогда все значения P_i равны
    double P = 0.0;
    for (int i = 0; i < n_fractions; ++i) {
        if (beta.has(i)) {
            double rho_i = (beta[i] / alpha[i]) * rho;
            P += m_materials[i]->pressure_rT(rho_i, T);
        }
    }
    P /= beta.nonzero();

    // Перенести объемные доли
    if (options.alpha) {
        *options.alpha = alpha;
    }

    return {P, T};
}

PairPT Materials::find_PT(double rho, double e, const Fractions &beta,
                          const Options &options) const {

    return find_PT_ver2(rho, e, beta, options);
}

double Materials::sound_speed_PT(double P, double T,
        const Fractions &beta, const Options &options) const {

    auto v = volume_PT(P, T, beta, {.deriv=true, .rho0=options.rho0, .alpha=options.alpha});
    auto e = energy_PT(P, T, beta, {.deriv=true, .rho0=options.rho0, .alpha=options.alpha});
    double inv_D = inv_J(v, e);

    double c2 = inv_D * sqr(v) * (e.dT + P * v.dT);
    return std::sqrt(c2);
}

} // namespace zephyr::phys