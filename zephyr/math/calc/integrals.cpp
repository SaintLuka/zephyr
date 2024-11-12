#include <cmath>

#include <zephyr/math/calc/integrals.h>

namespace zephyr::math {

// Одномерный интеграл
double integrate_trapz(const std::function<double(double)>& f, double a, double b, int n) {
    n = std::max(1, n);

    double h = (b - a) / n;

    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        res += f(a + (i + 0.5) * h);
    }
    return h * res;
}

double integrate_simpson(const std::function<double(double)>& f, double a, double b, int n) {
    const double c1 = 1.0 / 6.0;
    const double c2 = 1.0 / 3.0;
    const double c3 = 2.0 / 3.0;

    n = std::max(1, n);

    double h = (b - a) / n;

    double sum1 = f(a) + f(b);
    double sum2 = 0.0;
    double sum3 = f(a + 0.5 * h);
    for (int i = 1; i < n; ++i) {
        sum2 += f(a + i * h);
        sum3 += f(a + (i + 0.5) * h);
    }

    return h * (c1 * sum1 + c2 * sum2 + c3 * sum3);
}


double integrate_gauss(const std::function<double(double)>& f, double a, double b, int n) {
    const double coeff = 0.5 / std::sqrt(3.0);

    n = std::max(1, n);

    double h = (b - a) / n;
    double d = coeff * h;

    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        double x = a + (i + 0.5) * h;
        res += f(x - d) + f(x + d);
    }
    return 0.5 * h * res;
}

double integrate_gauss5(const std::function<double(double)>& f, double a, double b, int n) {
    const double w0 = 128.0 / 225;
    const double w1 = (322.0 + 13.0 * std::sqrt(70.0)) / 900.0;
    const double w2 = (322.0 - 13.0 * std::sqrt(70.0)) / 900.0;
    const double c1 = std::sqrt(5.0 - 2.0 * std::sqrt(10.0 / 7.0)) / 3.0;
    const double c2 = std::sqrt(5.0 + 2.0 * std::sqrt(10.0 / 7.0)) / 3.0;

    n = std::max(1, n);

    double h = (b - a) / n;
    double d1 = c1 * h;
    double d2 = c2 * h;

    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        double x = a + (i + 0.5) * h;
        res += w0 * f(x) +
               w1 * (f(x - d1) + f(x + d1)) +
               w2 * (f(x - d2) + f(x + d2));
    }
    return h * res;
}

/// n - число отрезков
double integrate(const std::function<double(double)>& f,
        double a, double b, int n, IntMethod method) {

    switch (method) {
        case IntMethod::TRAPZ:
            return integrate_trapz(f, a, b, n);
        case IntMethod::SIMPS:
            return integrate_simpson(f, a, b, n);
        case IntMethod::GAUSS:
            return integrate_gauss(f, a, b, n);
        default:
            return integrate_gauss5(f, a, b, n);
    }
}



} // namespace zephyr::math