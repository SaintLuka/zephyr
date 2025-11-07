// Протестируем параметрические кривые и сплайны, определенные в geom/curves.

#include <zephyr/geom/curves/bezier.h>
#include <zephyr/geom/curves/linear_spline.h>
#include <zephyr/geom/curves/cubic_spline.h>
#include <zephyr/geom/curves/lagrange.h>
#include <zephyr/utils/numpy.h>

#include <zephyr/utils/pyplot.h>

using namespace zephyr;
using namespace zephyr::geom;
using namespace zephyr::geom::curves;


// Уравнение спирали в полярных координатах r(phi)
Vector3d spiral(double angle) {
    const double phi = 0.5 * (1.0 + std::sqrt(5.0));
    double r = std::pow(phi, 2.0 * angle / M_PI);
    return {r * std::cos(angle), r * std::sin(angle), 0.0};
}

// Параметрическое задание эпитрохоиды
Vector3d epitr(double t) {
    double m = 0.2;
    double h = 0.3;
    double R = 1.0;
    double x = R * (m + 1) * std::cos(m * t) - h * std::cos((m + 1) * t);
    double y = R * (m + 1) * std::sin(m * t) - h * std::sin((m + 1) * t);
    return {x, y, 0.0};
}

// F -- параметрическая кривая
// t1, t2 -- область изменения параметра
// M -- число узлов интерполяции
// bezier -- построить кривую Безье
// chebyshev -- узлы Чебышёва для полинома Лагранжа
// periodic -- замкнутая периодическая кривая?
template <typename F>
void test(F&& f, double t1, double t2, int M, bool bezier, bool chebyshev, bool periodic) {
    // High-res curve
    int N = 2000;
    std::vector<double> phi_n = np::linspace(t1, t2, N);
    std::vector<Vector3d> vs_n(N);
    for (int i = 0; i < N; ++i) {
        vs_n[i] = f(phi_n[i]);
    }

    // Basic nodes
    std::vector<double> phi_m = np::linspace(t1, t2, M);
    std::vector<Vector3d> vs_m(M);
    for (int i = 0; i < M; ++i) {
        vs_m[i] = f(phi_m[i]);
    }

    utils::pyplot plt;

    plt.figure({.figsize={12.0, 8.0}});

    plt.grid(true);
    plt.set_aspect_equal();

    plt.plot(np::get_x(vs_n), np::get_y(vs_n), {.linestyle="dashed", .color="black", .label="Original"});

    PLagrange lg(vs_m, SplineBound::Crop, SplineBound::Crop,
                 chebyshev ? Parametrization::Chebyshev : Parametrization::Uniform);
    plt.plot(lg.xs(N), lg.ys(N), {.color="blue", .label="Lagrange"});

    PCubicSpline cs(vs_m, SplineBound::Crop,
                    periodic ? SplineBound::Periodic : SplineBound::Crop,
                    Parametrization::Uniform);
    plt.plot(cs.xs(N), cs.ys(N), {.color="green", .label="Cubic Spline"});

    if (bezier) {
        Bezier bz(vs_m);
        plt.plot(bz.xs(), bz.ys(), {.linestyle="dotted", .color="black", .label="Bezier basis"});
        plt.plot(bz.xs(N, 0.0, 1.05), bz.ys(N, 0.0, 1.05), {.color="orange", .label="Bezier curve"});
    }

    plt.plot(np::get_x(vs_m), np::get_y(vs_m), "k.");
    plt.legend();
    plt.tight_layout();
    plt.show();
}

int main() {
    test(spiral, 0.0, 3.0 * M_PI, 11, true, false, false);
    test(epitr, 0.0, 10.0 * M_PI, 25, false, true, true);

    return 0;
}