#include <iostream>

#include <zephyr/math/funcs.h>
#include <zephyr/math/calc/roots.h>
#include <zephyr/geom/sections.h>
#include <zephyr/math/calc/derivatives.h>

namespace zephyr::geom {

using namespace zephyr::math;

double quad_volume_fraction(const Vector3d &n, double p, double a, double b) {
    auto [xi, eta] = sorted(a * abs(n.x()), b * abs(n.y()));

    if (std::abs(p) <= 0.5 * (eta - xi)) {
        return 0.5 + p / eta;
    } else if (std::abs(p) < 0.5 * (eta + xi)) {
        return heav(p) - 0.5 * sign(p) * std::pow(std::abs(p) - 0.5 * (xi + eta), 2) / (xi * eta);
    } else {
        return heav(p);
    }
}

double quad_find_section(const Vector3d &n, double alpha, double a, double b) {
    auto [xi, eta] = sorted(a * abs(n.x()), b * abs(n.y()));

    if (alpha <= 0.5 * xi / eta) {
        return -0.5 * (xi + eta) + std::sqrt(2.0 * alpha * xi * eta);
    }
    else if (alpha < 1.0 - 0.5 * xi / eta) {
        return (alpha - 0.5) * eta;
    }
    else {
        return +0.5 * (xi + eta) - std::sqrt(2.0 * (1.0 - alpha) * xi * eta);
    }
}

namespace {

// Минимум и максимум, но не выходя за пределы [0, 1]
std::tuple<double, double> minmax_unit(double a1, double a2) {
    if (a1 < a2) {
        return {std::max(0.0, a1), std::min(a2, 1.0)};
    } else {
        return {std::max(0.0, a2), std::min(a1, 1.0)};
    }
}

}

double face_fraction_v3(double a1, double a2) {
    auto [a_min, a_max] = minmax_unit(a1, a2);

    if (a_min == 0.0) return 0.0;
    if (a_max == 1.0) return 1.0;

    if (a1 + a2 < 1.0) return a_max;
    if (a1 + a2 > 1.0) return a_min;

    return 0.5 * (a1 + a2);
}

double face_fraction_v5(double a1, double a2) {
    auto [a_min, a_max] = minmax_unit(a1, a2);

    if (a_min == 0.0) return 0.0;
    if (a_max == 1.0) return 1.0;

    return 0.5 * (a1 + a2);
}

double face_fraction_s(double a1, double a2) {
    auto [a_min, a_max] = minmax_unit(a1, a2);

    double a_sig = a_min / (1.0 - (a_max - a_min));

    // Случай a_min = 0, a_max = 1
    if (std::isnan(a_sig)) {
        a_sig = 0.5;
    }

    return  between(a_sig, a_min, a_max);
}

double face_fraction(double a1, double a2) {
    auto [a_min, a_max] = minmax_unit(a1, a2);

    if (3.0 * a_min >= a_max && a_min >= 3.0 * a_max - 2.0) {
        return 0.5 * (a_min + a_max);
    }
    else if (a_max <= 0.5 || a_min * (1.0 - a_max) >= std::pow(a_max - 0.5, 2)) {
        return -2.0 * a_min + 2.0 * std::sqrt(a_min * (a_min + a_max));
    }
    else if (a_min >= 0.5 || a_min * (1.0 - a_max) >= std::pow(a_min - 0.5, 2)) {
        return 3.0 - 2.0 * a_max - 2.0 * std::sqrt((1.0 - a_max) * (2.0 - a_min - a_max));
    }
    else {
        if (a_min == 0.0 && a_max == 1.0) {
            return 0.5;
        }
        return std::sqrt(a_min) / (std::sqrt(1.0 - a_max) + std::sqrt(a_min));
    }
}

double face_fraction_cos(double a1, double a2) {
    auto [a_min, a_max] = minmax_unit(a1, a2);

    double tg; // модуль тангенса
    if (3.0 * a_min >= a_max && a_min >= 3.0 * a_max - 2.0) {
        tg = 1.0 / std::abs(a1 - a2);
    } else if (a_max <= 0.5 || a_min * (1.0 - a_max) > std::pow(a_max - 0.5, 2)) {
        double ctg = 2.0 * a_max + 4.0 * a_min - 4.0 * std::sqrt(a_min * (a_min + a_max));
        tg = 1.0 / std::abs(ctg);
    } else if (a_min >= 0.5 || a_min * (1.0 - a_max) > std::pow(a_min - 0.5, 2)) {
        double ctg = 6.0 - 2.0 * a_min - 4.0 * a_max - 4.0 * std::sqrt((1.0 - a_max) * (2.0 - a_min - a_max));
        tg = 1.0 / std::abs(ctg);
    } else {
        tg = 2.0 * std::pow(std::sqrt(1.0 - a_max) + std::sqrt(a_min), 2);
    }

    return sign(a1 - a2) / std::sqrt(1.0 + tg * tg);
}

double average_flux(double alpha, double cos, double CFL) {
    auto [xi1, eta1] = sorted(std::abs(cos), std::sqrt(1.0 - cos * cos));
    auto [xi2, eta2] = sorted(CFL * std::abs(cos), std::sqrt(1.0 - cos * cos));

    // Фактически quad_find_section
    double P1;
    if (std::abs(2 * alpha - 1) <= 1.0 - xi1 / eta1) {
        P1 = (alpha - 0.5) * eta1;
    } else {
        P1 = sign(alpha - 0.5) * (0.5 * (xi1 + eta1) -
                                  std::sqrt((1.0 - std::abs(2.0 * alpha - 1.0)) * xi1 * eta1));
    }

    double P2 = P1 + 0.5 * (CFL - 1.0) * cos;

    // Фактически quad_volume_fraction
    if (std::abs(P2) <= 0.5 * (eta2 - xi2)) {
        return 0.5 + P2 / eta2;
    } else if (std::abs(P2) < 0.5 * (eta2 + xi2)) {
        return heav(P2) - 0.5 * sign(P2) * std::pow(std::abs(P2) - 0.5 * (xi2 + eta2), 2) / (xi2 * eta2);
    } else {
        return heav(P2);
    }
}

namespace {

// Положительная, монотонно возрастающая, кусочно-полиномиальная функция на отрезке.
// Параметры упорядочены 0 <= xi <= eta <= chi, одновременно не обращаются в ноль.
// Аргумент x ∈ [0, s], где s = (xi + eta + chi) / 2.
// Значения функции f(x) ∈ [0, 1].
double f_func(double x, double xi, double eta, double chi) {
    double s = 0.5 * (xi + eta + chi);

    if (x <= 0.5 * std::abs(xi + eta - chi)) {
        if (xi + eta <= chi) {
            // case 4.1. Квадратное сечение в центре
            return 2.0 * x / chi;
        } else {
            // case 4.2. Шестиугольное сечение в центре
            return x * (xi * eta + xi * chi + eta * chi - 0.5 * (pow(xi, 2) + pow(eta, 2) + pow(chi, 2)) - (2.0/ 3.0) * pow(x, 2)) / (xi * eta * chi);
        }
    }
    if (x <= 0.5 * (xi - eta + chi)) {
        // case 3:      0.5 |xi + eta - chi| < |x| <= 0.5 (xi - eta + chi)
        return -std::pow(x + s - chi, 3) / (3.0 * xi * eta * chi) + 2 * x / chi;
    }
    if (x <= 0.5 * (eta + chi - xi)) {
        // case 2:      0.5 (xi - eta + chi) < |x| <= 0.5 * (eta + chi - xi)
        return 1.0 + (3.0 * (x - s) * (-x + s - xi) - xi * xi) / (3.0 * eta * chi);
    }
    if (x <= s) {
        // case 1:      0.5 * (eta + chi - xi) < |x| <= s
        return 1.0 + pow(x - s, 3) / (3.0 * xi * eta * chi);
    }
    // case 0: |x| > s
    return 1.0;
}

// Тригонометрическая формула Виета для k-го корня кубического
// уравнения: y³ + py + q = 0
// Работает для случая трёх действительных корней
double viete(double p, double q, int k) {
    constexpr double c1 = 1.0 / 3.0;
    constexpr double c2 = 2.0 * M_PI / 3.0;

    double r = 2.0 * std::sqrt(-p / 3.0);
    double x = between(3.0 * q / (r * p), -1.0, 1.0);
    double phi = std::acos(x);
    return r * std::cos(c1 * phi + c2 * k);
}

// Обратная функция к f_func, также положительная и монотонно возрастающая на отрезке.
// Параметры упорядочены 0 <= xi <= eta <= chi, одновременно не обращаются в ноль.
// Аргумент y ∈ [0, 1], значения функции f(x) ∈ [0, s], где s = (xi + eta + chi) / 2.
double g_func(double y, double xi, double eta, double chi) {
    double s = 0.5 * (xi + eta + chi);

    if (xi + eta <= chi) {
        // case 4.1: квадратное сечение в центре
        if (y <= (chi - xi - eta) / chi) {
            return 0.5 * y * chi;
        }
    }
    else {
        // case 4.2: шестиугольное сечение в центре
        if (y <= (xi + eta - chi) / chi * (1.0 - std::pow(xi + eta - chi, 2) / (3.0 * xi * eta))) {
            double p = 0.75 * std::pow(xi + eta - chi, 2) - 3.0 * xi * eta;
            double q = 1.5 * xi * eta * chi * y;

            // Я проверил эти свойства в maple, но мало ли
            z_assert(p < 1.0e-12, "case 4.2: неверное предположение " + std::to_string(p));
            z_assert(pow(p / 3, 3) + pow(q / 2, 2) <= 0.0, "case 4.2: не случай трёх корней");

            // Корень должен оказаться на отрезке [0, 0.5 (xi + eta - chi)]
            z_assert(viete(p, q, 2) > -1.0e-12, "case 4.2: неверная оценка #1");
            z_assert(viete(p, q, 2) < 0.5 * (xi + eta - chi) + 1.0e-12, "case 4.2: неверная оценка #2");

            return viete(p, q, 2);
        }
    }
    if (y <= (xi - eta + chi) / chi - xi * xi / (3.0 * eta * chi)) {
        // case 3
        double p = -6.0 * xi * eta;
        double q = 3.0 * xi * eta * (xi + eta - chi * (1.0 - y));

        // Я проверил эти свойства в maple, но мало ли
        z_assert(q > -1.0e-12, "case 3: неверное предположение " + std::to_string(q));
        z_assert(pow(p / 3, 3) + pow(q / 2, 2) <= 0.0, "case 3: не случай трёх корней");

        // Корень должен оказаться на отрезке [0, xi]
        z_assert(viete(p, q, 2) > -1.0e-12, "case 3: неверная оценка #1");
        z_assert(viete(p, q, 2) < xi + 1.0e-12, "case 3: неверная оценка #2");

        return 0.5 * (chi - xi - eta) + viete(p, q, 2);
    }
    if (y <= 1.0 - xi * xi / (3.0 * eta * chi)) {
        // case 2
        return 0.5 * (chi + eta) - std::sqrt(chi * eta * (1.0 - y) - xi * xi / 12.0);
    }
    if (y <= 1.0) {
        // case 1
        return s - std::cbrt(3.0 * xi * eta * chi * (1.0 - y));
    }
    // case 0:
    return s;
}

}

double cube_volume_fraction(const Vector3d& n, double p, double a, double b, double c) {
    auto [xi, eta, chi] = sorted(std::abs(a * n.x()), std::abs(b * n.y()), std::abs(c * n.z()));
    return 0.5 + 0.5 * sign(p) * f_func(std::abs(p), xi, eta, chi);
}

double cube_find_section(const Vector3d& n, double alpha, double a, double b, double c) {
    auto [xi, eta, chi] = sorted(std::abs(a * n.x()), std::abs(b * n.y()), std::abs(c * n.z()));
    return sign(2.0 * alpha - 1.0) * g_func(std::abs(2.0 * alpha - 1.0), xi, eta, chi);
}

} // namespace zephyr::geom