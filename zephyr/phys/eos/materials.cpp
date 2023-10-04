#include <iostream>
#include <zephyr/phys/eos/materials.h>

namespace zephyr { namespace phys {

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

Eos& Materials::operator[](int idx){
    return *m_materials[idx];
}

const Eos& Materials::operator[](int idx) const {
    return *m_materials[idx];
}

dRdE Materials::pressure_re(double rho, double eps,
                            const Fractions& beta, const Options& options) const {
    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->pressure_re(rho, eps, options);
    }

    // Решаем v = sum_i beta_i v_i(P, T)
    double vol = 1.0 / rho;
    double P = std::isnan(options.T0) ? 1.0e5 : options.P0;
    double T = std::isnan(options.T0) ? 300.0 : options.T0;

    double err = 1.0;
    int counter = 0;
    while (err > 1.0e-10 && counter < 10) {
        auto v = volume_pt(P, T, beta);
        auto e = energy_pt(P, T, beta);

        // 1.0 / Якобиан
        double D = 1.0 / (v.dP * e.dT - v.dT * e.dP);

        double dP = D * ((vol - v.val) * e.dT - (eps - e.val) * v.dT);
        double dT = D * ((eps - e.val) * v.dP - (vol - v.val) * e.dP);

        err = std::abs(dP / P) + std::abs(dT / T);
        P += dP;
        T += dT;
        ++counter;
    }
    return P;
}

double Materials::energy_rp(double rho, double P,
        const Fractions& beta, const Options& options) const {
    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->energy_rp(rho, P, options);
    }

    double T = temperature_rp(rho, P, beta, options);
    return energy_pt(P, T, beta);
}

double Materials::sound_speed_re(double rho, double eps,
        const Fractions& beta, const Options& options) const {
    return NAN;
}

double Materials::sound_speed_rp(double rho, double P,
        const Fractions& beta, const Options& options) const {
    return NAN;
}

double Materials::pressure_rt(double rho, double T,
        const Fractions& beta, const Options& options) const {
    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->pressure_rt(rho, T, options);
    }

    // Решаем v = sum_i beta_i v_i(P, T)
    double vol = 1.0 / rho;
    double P = std::isnan(options.P0) ? 1.0e5 : options.P0;

    double err = 1.0;
    int counter = 0;
    while (err > 1.0e-10 && counter < 10) {
        auto v = volume_pt(P, T, beta);
        double dP = (vol - v.val) / v.dP;
        err = std::abs(dP / P);
        P += dP;
        ++counter;
    }
    return P;
}

double Materials::temperature_rp(double rho, double P,
        const Fractions &beta, const Options& options) const {
    // Случай одного материала
    int idx = beta.index();
    if (idx >= 0) {
        return m_materials[idx]->temperature_rp(rho, P, options);
    }

    // Решаем v = sum_i beta_i v_i(P, T)
    double vol = 1.0 / rho;
    double T = std::isnan(options.T0) ? 300.0 : options.T0;

    double err = 1.0;
    int counter = 0;
    while (err > 1.0e-10 && counter < 10) {
        auto v = volume_pt(P, T, beta);
        double dT = (vol - v.val) / v.dT;
        err = std::abs(dT / P);
        T += dT;
        ++counter;
    }
    return T;
}

dPdT Materials::volume_pt(double P, double T, const Fractions& beta) const {
    dPdT v = {0.0, 0.0, 0.0};
    for (int i = 0; i < size(); ++i) {
        if (beta.has(i)) {
            auto v_i = m_materials[i]->volume_pt(P, T, {.deriv = true});
            v.val += beta[i] * v_i.val;
            v.dP  += beta[i] * v_i.dP;
            v.dT  += beta[i] * v_i.dT;
        }
    }
    return v;
}

dPdT Materials::energy_pt(double P, double T, const Fractions& beta) const {
    dPdT e = {0.0, 0.0, 0.0};
    for (int i = 0; i < size(); ++i) {
        if (beta.has(i)) {
            auto e_i = m_materials[i]->energy_pt(P, T, {.deriv = true});
            e.val += beta[i] * e_i.val;
            e.dP  += beta[i] * e_i.dP;
            e.dT  += beta[i] * e_i.dT;
        }
    }
    return e;
}

StiffenedGas Materials::stiffened_gas(double rho, double P,
                   const Fractions& beta, const Options& options) const {
    return StiffenedGas(NAN, NAN, NAN, NAN);
}

}}