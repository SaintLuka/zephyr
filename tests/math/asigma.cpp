#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <zephyr/math/funcs.h>
#include <zephyr/geom/vector.h>

#include <zephyr/utils/matplotlib.h>
namespace plt = zephyr::utils::matplotlib;

using namespace zephyr;

void test_heaviside() {
    auto x = geom::linspace(-1.0, 3.0, 501);

    std::vector<double> H_s(x.size());
    std::vector<double> H_p2(x.size());
    std::vector<double> H_p4(x.size());
    std::vector<double> H_r3(x.size());
    std::vector<double> H_r6(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        H_s[i] = math::heav_s(x[i]);
        H_p2[i] = math::heav_p<2>(x[i]);
        H_p4[i] = math::heav_p(x[i]);
        H_r6[i] = math::heav_r<6>(x[i]);
        H_r3[i] = math::heav_r(x[i]);
    }

    plt::figure_size(9.0, 3.0, 200);

    plt::grid(true);

    plt::xlim(-1.0, 3.0);
    plt::title("Гладкая функция Хевисайда");

    plt::plot(x, H_s,  {{"label", "$\\eta_{sin}(x)$"}, {"linestyle", "dotted"}});
    plt::plot(x, H_p2, {{"label", "$\\eta_p(x), n=2$"}, {"linestyle", "dashed"}});
    plt::plot(x, H_p4, {{"label", "$\\eta_p(x), n=4$"}, {"linestyle", "dashed"}});
    plt::plot(x, H_r3, {{"label", "$\\eta_r(x), n=3$"}, {"linestyle", "solid"}});
    plt::plot(x, H_r6, {{"label", "$\\eta_r(x), n=6$"}, {"linestyle", "solid"}});
    plt::legend();

    plt::tight_layout();
    plt::show();
}

void test_sigmoid() {
    auto x = geom::linspace(-3.0, 3.0, 601);

    std::vector<double> sgn_s(x.size());
    std::vector<double> sgn_p2(x.size());
    std::vector<double> sgn_p4(x.size());
    std::vector<double> sgn_r3(x.size());
    std::vector<double> sgn_r6(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        sgn_s[i]  = math::sign_s(x[i]);
        sgn_p2[i] = math::sign_p<2>(x[i]);
        sgn_p4[i] = math::sign_p(x[i]);
        sgn_r6[i] = math::sign_r<6>(x[i]);
        sgn_r3[i] = math::sign_r(x[i]);
    }

    plt::figure_size(9.0, 3.0, 200);

    plt::grid(true);

    plt::xlim(-3.0, 3.0);
    plt::title("Гладкая функция знака / сигмоида");

    plt::plot(x, sgn_s,  {{"label", "sgn$_{sin}(x)$"}, {"linestyle", "dotted"}});
    plt::plot(x, sgn_p2, {{"label", "sgn$_p(x), n=2$"}, {"linestyle", "dashed"}});
    plt::plot(x, sgn_p4, {{"label", "sgn$_p(x), n=4$"}, {"linestyle", "dashed"}});
    plt::plot(x, sgn_r3, {{"label", "sgn$_r(x), n=3$"}, {"linestyle", "solid"}});
    plt::plot(x, sgn_r6, {{"label", "sgn$_r(x), n=6$"}, {"linestyle", "solid"}});
    plt::legend();

    plt::tight_layout();
    plt::show();
}

double a_sig_v3(double a1, double a2) {
    auto [a_min, a_max] = math::minmax(a1, a2);

    if (a_min == 0.0) return 0.0;
    if (a_max == 1.0) return 1.0;

    if (a1 + a2 < 1.0) return a_max;
    if (a1 + a2 > 1.0) return a_min;

    return 0.5 * (a1 + a2);
}

double a_sig_v5(double a1, double a2) {
    auto [a_min, a_max] = math::minmax(a1, a2);

    if (a_min == 0.0) return 0.0;
    if (a_max == 1.0) return 1.0;

    return 0.5 * (a1 + a2);
}

double a_sig_v3s(double a1, double a2) {
    auto [a_min, a_max] = math::minmax(a1, a2);

    if (a_min == 0.0) return 0.0;
    if (a_max == 1.0) return 1.0;

    double gamma = a1 + a2 - 1.0;
    double delta = a_max - a_min;

    double H = math::heav_p(1.0 - delta - std::abs(gamma), 0.3 * delta);
    double S = math::sign_p(gamma, 0.2 * delta);

    double theta = 1.0 + H * (std::abs(gamma) - delta - 1.0);
    return 0.5 * (1.0 + S * theta);
}

double a_sig_v5s(double a1, double a2) {
    auto [a_min, a_max] = math::minmax(a1, a2);

    if (a_min == 0.0) return 0.0;
    if (a_max == 1.0) return 1.0;

    double gamma = a1 + a2 - 1.0;
    double delta = a_max - a_min;

    double H = math::heav_p(1.0 - delta - std::abs(gamma), 0.3 * delta);

    double theta = 1.0 + H * (std::abs(gamma) - 1.0);
    return 0.5 * (1.0 + math::sign(gamma) * theta);
}

int main() {
    //test_heaviside();
    //test_sigmoid();

    int n = 101;

    auto[a1, a2] = geom::meshgrid(
            geom::linspace(0.0, 1.0, n),
            geom::linspace(0.0, 1.0, n));

    std::vector<std::vector<double>> f1(n, std::vector<double>(n));
    std::vector<std::vector<double>> f2(n, std::vector<double>(n));
    std::vector<std::vector<double>> f3(n, std::vector<double>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            f1[i][j] = a_sig_v3(a1[i][j], a2[i][j]);
            f2[i][j] = a_sig_v3s(a1[i][j], a2[i][j]);

            double v1 = a_sig_v5s(a1[i][j], a1[i][j]);
            double v2 = a_sig_v5s(1.0 - a1[i][j], 1.0 - a2[i][j]);
            f3[i][j] = f1[i][j] - f2[i][j];
        }
    }

    plt::plot_surface(a1, a2, f1, {{"cmap", "jet"}});
    plt::plot_surface(a1, a2, f2, {{"cmap", "jet"}});
    plt::plot_surface(a1, a2, f3, {{"cmap", "jet"}});
    plt::show();

    return 0;
}