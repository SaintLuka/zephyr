#include <vector>

#include <zephyr/phys/tests/sedov.h>

namespace zephyr::phys {

Sedov3D::Sedov3D(double _gamma, double _rho0, double _E) {
    gamma = _gamma;
    rho0  = _rho0;
    E     = _E;

    // Константы
    double kek = (13 * gamma * gamma - 7 * gamma + 12);
    double nu1 = -kek / ((2 * gamma + 1) * (3 * gamma - 1));
    double nu2 = kek / (5 * (2 - gamma) * (3 * gamma - 1));
    double nu3 = kek / (5 * (2 - gamma) * (gamma - 1) * (3 * gamma - 1));
    double nu4 = (2 * gamma + 1) / (gamma - 1);
    double nu5 = -2 / (2 - gamma);
    double nu6 = 3 / (gamma - 1);
    double nu7 = 6 / (5 * (gamma - 1));

    double V_min = 1.0 / gamma;
    double V_max = 2.0 / (gamma + 1.0);

    // Строим линейный сплайн. Здесь важно правильно выбрать
    // узлы интерполяции, поскольку в нуле сингулярная точка.
    // Узлы выбираются так, что формируется почти равномерный
    // массив xi[i], получено так V(xi) ~ V_min + C s^nu4.
    int n_nodes = 1000;
    m_xi.resize(n_nodes);
    m_Vs.resize(n_nodes);
    m_Zs.resize(n_nodes);
    m_Gs.resize(n_nodes);
    m_Ps.resize(n_nodes);

    std::vector<double> subint(n_nodes);
    for (int i = 0; i < n_nodes; ++i) {
        double s = double(i) / (n_nodes - 1);
        double q = std::pow(s, nu4);

        double V = V_min + (V_max - V_min) * q;

        double F1 = V / V_max;
        double F2 = std::pow((gamma - 1.0) / (gamma + 1.0), -0.5 * gamma) * (1 - V);
        double F3 = (gamma + 1.0) / (7.0 - gamma) * (5.0 - (3.0 * gamma - 1.0) * V);

        double xi = std::pow(std::pow(F3, nu1) / std::pow(F1, 2), 0.2) * s;

        double Z = std::pow(V, 2) * (1 - V) / (V_min * V_max * q);

        double G = std::pow(F1, nu7) *
                   std::pow(F2, nu5) *
                   std::pow(F3, nu3) *
                   std::pow(xi, nu6);

        double P = std::pow(xi, 2) * Z * G;

        m_xi[i] = xi;
        m_Vs[i] = V;
        m_Zs[i] = Z;
        m_Gs[i] = G;
        m_Ps[i] = P;

        subint[i] = std::pow(V, 2) *
                    std::pow(F1, 0.2) *
                    std::pow(F2, nu5) *
                    std::pow(F3, nu2) *
                    std::pow(xi, 2);
    }

    double b_int = 0.0;
    for (int i = 1; i < n_nodes; ++i) {
        b_int += 0.5 * (subint[i] + subint[i - 1]) * (m_xi[i] - m_xi[i - 1]);
    }

    double K = 16.0 * M_PI / 25.0;
    beta = std::pow(K * b_int, -0.2);
}

// Линейная интерполяция
double interpolation(const std::vector<double>& xs, const std::vector<double>& ys, double x) {
    x = std::max(0.0, std::min(x, 1.0));
    int n = std::distance(xs.begin(), std::lower_bound(xs.begin(), xs.end(), x));
    n = std::max(0, n - 1);
    double x1 = xs[n];
    double x2 = xs[n + 1];
    double y1 = ys[n];
    double y2 = ys[n + 1];
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

double Sedov3D::V_xi(double xi) const {
    return interpolation(m_xi, m_Vs, xi);
}

double Sedov3D::Z_xi(double xi) const {
    return interpolation(m_xi, m_Zs, xi);
}

double Sedov3D::G_xi(double xi) const {
    return interpolation(m_xi, m_Gs, xi);
}

double Sedov3D::P_xi(double xi) const {
    return interpolation(m_xi, m_Ps, xi);
}

double Sedov3D::r_shock(double t) const {
    return beta * std::pow(E * t * t / rho0, 0.2);
}

double Sedov3D::time_by_radius(double r) const {
    return std::sqrt(rho0 / E * std::pow(r / beta, 5));
}

double Sedov3D::c_squared(double r, double t) const {
    double xi = r / r_shock(t);
    return xi <= 1.0 ? std::pow(0.4 * r / t, 2) * Z_xi(xi) : 0.0;
}

double Sedov3D::density(double r, double t) const {
    double xi = r / r_shock(t);
    return xi <= 1.0 ? rho0 * G_xi(xi) : rho0;
}

double Sedov3D::velocity(double r, double t) const {
    double xi = r / r_shock(t);
    return xi <= 1.0 ? (0.4 * r / t) * V_xi(xi) : 0.0;
}

double Sedov3D::pressure(double r, double t) const {
    double xi = r / r_shock(t);
    return xi <= 1.0 ? rho0 * std::pow(0.4 * r_shock(t) / t, 2) * P_xi(xi) / gamma : 0.0;
}

double Sedov3D::energy(double r, double t) const {
    return c_squared(r, t) / (gamma * (gamma - 1.0));
}

} // namespace zephyr::phys