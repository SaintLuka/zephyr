#include <zephyr/utils/numpy.h>
#include <zephyr/utils/pyplot.h>
#include <zephyr/phys/tests/sedov.h>

using namespace zephyr;
using phys::Sedov3D;

int main() {
    Sedov3D test({.gamma=1.4, .rho0=1.0, .E=1.0});

    double t0 = test.time_by_radius(0.6);

    int n = 300;
    auto r = np::linspace(0.0, 1.0, n);

    auto rho = np::zeros(n);
    auto P   = np::zeros(n);
    auto v   = np::zeros(n);
    auto e   = np::zeros(n);

    for (int i = 0; i < n; ++i) {
        rho[i] = test.density (r[i], t0);
        P[i]   = test.pressure(r[i], t0);
        v[i]   = test.velocity(r[i], t0);
        e[i]   = test.energy  (r[i], t0);
    }

    utils::pyplot plt;

    plt.figure({.figsize={9.0, 6.0}, .dpi=150});

    plt.subplot(2, 2, 1);
    plt.title("$\\rho(r, t_0)$");
    plt.grid(true);
    plt.xlim(0.0, 1.0);
    plt.plot(r, rho, "k");

    plt.subplot(2, 2, 2);
    plt.title("$P(r, t_0)$");
    plt.grid(true);
    plt.xlim(0.0, 1.0);
    plt.plot(r, P, "k");

    plt.subplot(2, 2, 3);
    plt.title("$v(r, t_0)$");
    plt.grid(true);
    plt.xlim(0.0, 1.0);
    plt.plot(r, v, "k");

    plt.subplot(2, 2, 4);
    plt.title("$e(r, t_0)$");
    plt.grid(true);
    plt.xlim(0.0, 1.0);
    plt.semilogy(r, e, "k");

    plt.tight_layout();
    plt.show();

    return 0;
}