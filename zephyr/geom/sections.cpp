#include <iostream>

#include <zephyr/math/funcs.h>
#include <zephyr/geom/sections.h>

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

double cube_volume_fraction(const Vector3d& n, double p, double a, double b, double c) {
    auto [xi, eta, chi] = sorted(std::abs(a * n.x()), std::abs(b * n.y()), std::abs(c * n.z()));
    double s = 0.5 * (xi + eta + chi);

    //if (xi + eta < chi) { std::cout << "quad case\n"; }
    //else { std::cout << "hex case\n"; }

    double f = NAN;
    double ap = std::abs(p);

    if (ap <= 0.5 * std::abs(xi + eta - chi)) {
        if (xi + eta < chi) {
            // case 4.1. Квадратное сечение в центре
            f = 2.0 * ap / chi;
        } else {
            // case 4.2. Шестиугольное сечение в центре
            f = ap * (xi * eta + xi * chi + eta * chi - 0.5 * (pow(xi, 2) + pow(eta, 2) + pow(chi, 2)) - (2.0/ 3.0) * pow(ap, 2)) / (xi * eta * chi);
        }
    }
    else if (ap <= 0.5 * (xi - eta + chi)) {
        // case 3:      0.5 |xi + eta - chi| < |p| <= 0.5 (xi - eta + chi)
        f = 1.0 + (-pow(s - ap, 3) + pow(s - ap - xi, 3) + pow(s - ap - eta, 3)) / (3.0 * xi * eta * chi);
    }
    else if (ap <= 0.5 * (eta + chi - xi)) {
        // case 2:      0.5 (xi - eta + chi) < |p| <= 0.5 * (eta + chi - xi)
        f = 1.0 + (3.0 * (ap - s) * (-ap + s - xi) - xi * xi) / (3.0 * eta * chi);
    }
    else if (ap <= s) {
        // case 1:      0.5 * (eta + chi - xi) < |p| <= s
        f = 1.0 + pow(ap - s, 3) / (3.0 * xi * eta * chi);
    }
    else {
        // case 0:      |p| > s
        f = 1.0;
    }
    return 0.5 + 0.5 * sign(p) * f;
}

double cube_find_section(const Vector3d& n, double alpha, double a, double b, double c) {
    return NAN;
}

} // namespace zephyr::geom