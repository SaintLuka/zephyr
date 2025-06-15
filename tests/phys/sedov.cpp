#include <vector>

#include <zephyr/utils/matplotlib.h>
#include <zephyr/phys/tests/sedov.h>

namespace plt = zephyr::utils::matplotlib;

using zephyr::phys::SedovBlast3D;

int main() {
    SedovBlast3D test({.gamma=1.4, .rho0=1.0, .E=1.0});

    double t0 = test.time_by_radius(0.6);

    int n = 300;
    std::vector<double> r(n);

    std::vector<double> rho(n);
    std::vector<double> P(n);
    std::vector<double> v(n);
    std::vector<double> e(n);

    for (int i = 0; i < n; ++i) {
        r[i] = double(i) / (n - 1);

        rho[i] = test.density (r[i], t0);
        P[i]   = test.pressure(r[i], t0);
        v[i]   = test.velocity(r[i], t0);
        e[i]   = test.energy  (r[i], t0);
    }

    plt::figure_size(9.0, 6.0, 170);

    plt::subplot(2, 2, 1);
    plt::title("$\\rho(r, t_0)$");
    plt::grid(true);
    plt::xlim(0.0, 1.0);
    plt::plot(r, rho, "k");

    plt::subplot(2, 2, 2);
    plt::title("$P(r, t_0)$");
    plt::grid(true);
    plt::xlim(0.0, 1.0);
    plt::plot(r, P, "k");

    plt::subplot(2, 2, 3);
    plt::title("$v(r, t_0)$");
    plt::grid(true);
    plt::xlim(0.0, 1.0);
    plt::plot(r, v, "k");

    plt::subplot(2, 2, 4);
    plt::title("$e(r, t_0)$");
    plt::grid(true);
    plt::xlim(0.0, 1.0);
    plt::semilogy(r, e, "k");

    plt::tight_layout();
    plt::show();

    return 0;
}