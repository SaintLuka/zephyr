#include <iostream>

#include <zephyr/geom/sections.h>
#include <zephyr/math/funcs.h>
#include <zephyr/math/calc/roots.h>
#include <zephyr/math/calc/derivatives.h>
#include <zephyr/mesh/side.h>

namespace zephyr::geom {

using namespace zephyr::math;

double quad_volume_fraction(double p, const Vector3d &n, double A, double B) {
    auto [a, b] = sorted(A * abs(n.x()), B * abs(n.y()));

    if (std::abs(p) <= 0.5 * (b - a)) {
        return 0.5 + p / b;
    } else if (std::abs(p) < 0.5 * (b + a)) {
        return heav(p) - 0.5 * sign(p) * std::pow(std::abs(p) - 0.5 * (a + b), 2) / (a * b);
    } else {
        return heav(p);
    }
}

double quad_find_section(double alpha, const Vector3d &n, double A, double B) {
    auto [a, b] = sorted(A * abs(n.x()), B * abs(n.y()));

    if (alpha <= 0.5 * a / b) {
        return -0.5 * (a + b) + std::sqrt(2.0 * alpha * a * b);
    }
    else if (alpha < 1.0 - 0.5 * a / b) {
        return (alpha - 0.5) * b;
    }
    else {
        return +0.5 * (a + b) - std::sqrt(2.0 * (1.0 - alpha) * a * b);
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
    auto [a1, b1] = sorted(std::abs(cos), std::sqrt(1.0 - cos * cos));
    auto [a2, b2] = sorted(CFL * std::abs(cos), std::sqrt(1.0 - cos * cos));

    // Фактически quad_find_section
    double P1;
    if (std::abs(2 * alpha - 1) <= 1.0 - a1 / b1) {
        P1 = (alpha - 0.5) * b1;
    } else {
        P1 = sign(alpha - 0.5) * (0.5 * (a1 + b1) -
                                  std::sqrt((1.0 - std::abs(2.0 * alpha - 1.0)) * a1 * b1));
    }

    double P2 = P1 + 0.5 * (CFL - 1.0) * cos;

    // Фактически quad_volume_fraction
    if (std::abs(P2) <= 0.5 * (b2 - a2)) {
        return 0.5 + P2 / b2;
    } else if (std::abs(P2) < 0.5 * (b2 + a2)) {
        return heav(P2) - 0.5 * sign(P2) * std::pow(std::abs(P2) - 0.5 * (a2 + b2), 2) / (a2 * b2);
    } else {
        return heav(P2);
    }
}

namespace {

// Положительная, монотонно возрастающая, кусочно-полиномиальная функция на отрезке.
// Параметры упорядочены 0 <= a <= b <= c, одновременно не обращаются в ноль.
// Аргумент x ∈ [0, s], где s = (a + b + c) / 2.
// Значения функции f(x) ∈ [0, 1].
double f_func(double x, double a, double b, double c) {
    double s = 0.5 * (a + b + c);

    if (x <= std::abs(s - c)) {
        if (a + b <= c) {
            // case 4.1. Квадратное сечение в центре
            return 2.0 * x / c;
        } else {
            // case 4.2. Шестиугольное сечение в центре
            return 2.0 * x * (a * b - std::pow(s - c, 2) - pow(x, 2) / 3.0) / (a * b * c);
        }
    }
    if (x < s - b) {
        // case 3:      0.5 |a + b - c| < |x| <= 0.5 (a - b + c)
        return 2 * x / c - std::pow(x + s - c, 3) / (3.0 * a * b * c);
    }
    if (x <= s - a) {
        // case 2:      0.5 (a - b + c) < |x| <= 0.5 * (b + c - a)
        return 1.0 - (3.0 * (s - x) * (s - x - a) + a * a) / (3.0 * b * c);
    }
    if (x <= s) {
        // case 1:      0.5 * (b + c - a) < |x| <= s
        return 1.0 - pow(s - x, 3) / (3.0 * a * b * c);
    }
    // case 0: |x| > s
    return 1.0;
}

vf_grad_t f_func_grad(double x, double a, double b, double c) {
    double s = 0.5 * (a + b + c);

    if (x <= std::abs(s - c)) {
        if (a + b <= c) {
            // case 4.1. Квадратное сечение в центре
            return {2.0 / c, Vector3d{0.0, 0.0, -2.0 * x / (c * c)}};
        } else {
            // case 4.2. Шестиугольное сечение в центре
            double dx = 2.0 / c - 2.0 * (pow(s - c, 2) + x * x) / (a * b * c);
            double XX = 2.0 * x / (a * b * c);
            double da = (x*x / 3.0 - (s - b) * (s - c)) * XX / a;
            double db = (x*x / 3.0 - (s - a) * (s - c)) * XX / b;
            double dc = (x*x / 3.0 - (s - a) * (s - b)) * XX / c;
            return {dx, Vector3d{da, db, dc}};
        }
    }
    if (x < s - b) {
        // case 3:      0.5 |a + b - c| < |x| <= 0.5 (a - b + c)
        double dx = 2.0 / c - pow(x + s - c, 2) / (a * b * c);
        double XX = pow(x + s - c, 2) / (3.0 * a * b * c);
        double da = (x - a + 0.5 * (b - c)) * XX / a;
        double db = (x - b + 0.5 * (a - c)) * XX / b;
        double dc = -2.0*x/(c*c) + (x + s + 0.5 * c) * XX / c;
        return {dx, Vector3d{da, db, dc}};
    }
    if (x <= s - a) {
        // case 2:      0.5 (a - b + c) < |x| <= 0.5 * (b + c - a)
        double dx = (b + c - 2 * x)/(b * c);
        double da = -a / (6.0 * b * c);
        double db = (a*a - 3*b*b + 3*c*c + 12*x*(x - c))/(12 * b * b * c);
        double dc = (a*a + 3*b*b - 3*c*c + 12*x*(x - b))/(12 * b * c * c);
        return {dx, Vector3d{da, db, dc}};
    }
    if (x <= s) {
        // case 1:      0.5 * (b + c - a) < |x| <= s
        double dx = pow(s - x, 2) / (a * b * c);
        double da = dx * ((s - x) / (3.0 * a) - 0.5);
        double db = dx * ((s - x) / (3.0 * b) - 0.5);
        double dc = dx * ((s - x) / (3.0 * c) - 0.5);
        return {dx, Vector3d{da, db, dc}};
    }
    // case 0: |x| > s
    return {0.0, Vector3d::Zero()};
}

int f_func_case(double x, double a, double b, double c) {
    double s = 0.5 * (a + b + c);

    if (x <= std::abs(s - c)) {
        if (a + b <= c) {
            // case 4.1. Квадратное сечение в центре
            return 5;
        } else {
            // case 4.2. Шестиугольное сечение в центре
            return 4;
        }
    }
    if (x < s - b) {
        // case 3:      0.5 |a + b - c| < |x| <= 0.5 (a - b + c)
        return 3;
    }
    if (x <= s - a) {
        // case 2:      0.5 (a - b + c) < |x| <= 0.5 * (b + c - a)
        return 2;
    }
    if (x <= s) {
        // case 1:      0.5 * (b + c - a) < |x| <= s
        return 1;
    }
    // case 0: |x| > s
    return 0;
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
// Параметры упорядочены 0 <= a <= b <= c, одновременно не обращаются в ноль.
// Аргумент y ∈ [0, 1], значения функции f(x) ∈ [0, s], где s = (a + b + c) / 2.
double g_func(double y, double a, double b, double c) {
    double s = 0.5 * (a + b + c);

    if (a + b <= c) {
        // case 4.1: квадратное сечение в центре
        if (y <= (c - a - b) / c) {
            return 0.5 * y * c;
        }
    }
    else {
        // case 4.2: шестиугольное сечение в центре
        if (y <= (a + b - c) / c * (1.0 - std::pow(a + b - c, 2) / (3.0 * a * b))) {
            double p = 3.0 * (std::pow(s - c, 2) - a * b);
            double q = 1.5 * a * b * c * y;

            // Я проверил эти свойства в maple, но мало ли
            z_assert(p < 1.0e-12, "case 4.2: неверное предположение " + std::to_string(p));
            z_assert(pow(p / 3, 3) + pow(q / 2, 2) <= 0.0, "case 4.2: не случай трёх корней");

            // Корень должен оказаться на отрезке [0, s - c]
            z_assert(viete(p, q, 2) > -1.0e-12, "case 4.2: неверная оценка #1");
            z_assert(viete(p, q, 2) < s - c + 1.0e-12, "case 4.2: неверная оценка #2");

            return viete(p, q, 2);
        }
    }
    if (y < (a - b + c) / c - a * a / (3.0 * b * c)) {
        // case 3
        double p = -6.0 * a * b;
        double q = 6.0 * a * b * (s - c * (1.0 - 0.5 * y));

        // Я проверил эти свойства в maple, но мало ли
        z_assert(q > -1.0e-12, "case 3: неверное предположение " + std::to_string(q));
        z_assert(pow(p / 3, 3) + pow(q / 2, 2) <= 0.0, "case 3: не случай трёх корней");

        // Корень должен оказаться на отрезке [0, a]
        z_assert(viete(p, q, 2) > -1.0e-12, "case 3: неверная оценка #1");
        z_assert(viete(p, q, 2) < a + 1.0e-12, "case 3: неверная оценка #2");

        return 0.5 * (c - a - b) + viete(p, q, 2);
    }
    if (y <= 1.0 - a * a / (3.0 * b * c)) {
        // case 2
        return 0.5 * (b + c) - std::sqrt(c * b * (1.0 - y) - a * a / 12.0);
    }
    if (y <= 1.0) {
        // case 1
        return s - std::cbrt(3.0 * a * b * c * (1.0 - y));
    }
    // case 0:
    return s;
}

}

double cube_volume_fraction(double p, const Vector3d& n, double A, double B, double C) {
    auto [a, b, c] = sorted(std::abs(A * n.x()), std::abs(B * n.y()), std::abs(C * n.z()));
    return 0.5 + 0.5 * sign(p) * f_func(std::abs(p), a, b, c);
}

vf_grad_t cube_volume_fraction_grad(double p, const Vector3d& n, double A, double B, double C) {
    double a = std::abs(A * n.x());
    double b = std::abs(B * n.y());
    double c = std::abs(C * n.z());

    std::array<bool, 3> swaps{false, false, false};

    // Сортируем по возрастанию, запоминаем перестановки
    if (a > b) { swaps[0] = true; std::swap(a, b); }
    if (b > c) { swaps[1] = true; std::swap(b, c); }
    if (a > b) { swaps[2] = true; std::swap(a, b); }

    auto grad = f_func_grad(std::abs(p), a, b, c);

    // Обратная перестановка
    if (swaps[2]) { std::swap(grad.dn[0], grad.dn[1]); }
    if (swaps[1]) { std::swap(grad.dn[1], grad.dn[2]); }
    if (swaps[0]) { std::swap(grad.dn[0], grad.dn[1]); }

    // Масштабирование
    grad.dp *= 0.5;

    double sign_p = p < 0.0 ? -0.5 : 0.5;
    grad.dn[0] *= sign_p * A * sign(n.x());
    grad.dn[1] *= sign_p * B * sign(n.y());
    grad.dn[2] *= sign_p * C * sign(n.z());

    return grad;

    // Сравнить с численным градиентом (в режиме debug)
    debug_code {
        if (n.norm() == 0.0) {
            throw std::invalid_argument("zero normal");
        }

        vf_grad_t grad2;

        auto vf_p  = [A, B, C, n, p](double t) -> double { return cube_volume_fraction(t, n, A, B, C); };
        auto vf_nx = [A, B, C, n, p](double t) -> double { return cube_volume_fraction(p, Vector3d{t, n.y(), n.z()}, A, B, C); };
        auto vf_ny = [A, B, C, n, p](double t) -> double { return cube_volume_fraction(p, Vector3d{n.x(), t, n.z()}, A, B, C); };
        auto vf_nz = [A, B, C, n, p](double t) -> double { return cube_volume_fraction(p, Vector3d{n.x(), n.y(), t}, A, B, C); };

        double h = 1.0e-5;
        grad2.dp    = derivative<1>(vf_p, p, h);
        grad2.dn[0] = derivative<1>(vf_nx, n.x(), h);
        grad2.dn[1] = derivative<1>(vf_ny, n.y(), h);
        grad2.dn[2] = derivative<1>(vf_nz, n.z(), h);

        vf_grad_t err = {
            .dp = std::abs(grad.dp - grad2.dp),
            .dn = (grad.dn - grad2.dn).cwiseAbs()
        };

        if (std::max(err.dp, err.dn.norm()) > 1.0e-5) {
            double s = 0.5 * (a + b + c);
            std::cout << "params: " << p << " " << n.transpose() << ";\n";
            std::cout << "a b c: " << a << " " << b << " " << c << "; " << p + s - b << "\n";
            std::cout << "vals: " << err.dp << " " << err.dn.transpose() << "\n";
            std::cout << "pairs:\n";
            std::cout << "case: " << f_func_case(std::abs(p), a, b, c) << "\n";
            throw std::runtime_error("volume_fraction_grad error");
        }
    }

    return grad;
}

double cube_find_section(double alpha, const Vector3d& n, double A, double B, double C) {
    auto [a, b, c] = sorted(std::abs(A * n.x()), std::abs(B * n.y()), std::abs(C * n.z()));
    return sign(2.0 * alpha - 1.0) * g_func(std::abs(2.0 * alpha - 1.0), a, b, c);
}

int cube_section_case(double p, const Vector3d& n, double A, double B, double C) {
    auto [a, b, c] = sorted(std::abs(A * n.x()), std::abs(B * n.y()), std::abs(C * n.z()));
    return (p > 0.0 ? 1 : -1) * f_func_case(std::abs(p), a, b, c);
}

inline double get_ae(double p, const Vector3d& n) {
    if (std::abs(n.z()) > 1.0e-15) {
        double z = (p - 0.5 * n.x() - 0.5 * n.y()) / n.z();
        return between(0.5 + z, 0.0, 1.0);
    }
    return n.x() + n.y() > 0.0 ? 0.0 : 1.0;
}

double edge_fraction(double a0, double a1, double a2) {
    constexpr double eps_a = 1.0e-12;
    if (a1 <= eps_a || a2 <= eps_a) {
        return 0.0;
    }
    if (a1 >= 1.0 - eps_a || a2 >= 1.0 - eps_a) {
        return 1.0;
    }
    if (a0 <= eps_a || a0 >= 1.0 - eps_a) {
        // return 0 или 1
        return 0.5 + 0.5 * math::sign_p(a1 + a2 - 1.0, eps_a);
    }

    // Остались случаи:
    //  eps_a < a0 < 1 - eps_a,
    //  eps_a < a1 < 1 - eps_a,
    //  eps_a < a2 < 1 - eps_a

    // Это как первая итерация в геометрическом методе, работает отлично!
    // Для простого случая сразу дает ответ.
    Vector3d n = {a0 - a1, a0 - a2, 1}; n.normalize();
    double p = (a0 - 0.5) * n.z();
    double ae = get_ae(p, n);

    constexpr int max_iters = 20;
    constexpr double eps_r = 1.0e-12;
    constexpr double eps_p = 1.0e-12;
    constexpr double eps_n = 1.0e-12;
    constexpr double eps_ae = 1.0e-12;

    for (int counter = 0; counter < max_iters; ++counter) {
        double A0 = cube_volume_fraction(p, n);
        double A1 = cube_volume_fraction(p - n.x(), n);
        double A2 = cube_volume_fraction(p - n.y(), n);

        Vector4d F = { a0 - A0, a1 - A1, a2 - A2, 0.0 };

        // Условие выхода по максимальной невязке
        double err_1 = F.cwiseAbs().maxCoeff();
        if (err_1 < eps_r) {
            break;
        }

        vf_grad_t D0 = cube_volume_fraction_grad(p, n);
        vf_grad_t D1 = cube_volume_fraction_grad(p - n.x(), n);
        vf_grad_t D2 = cube_volume_fraction_grad(p - n.y(), n);

        Matrix4d M;
        M << D0.dp, D0.dn.x(),         D0.dn.y(),         D0.dn.z(),
             D1.dp, D1.dn.x() - D1.dp, D1.dn.y(),         D1.dn.z(),
             D2.dp, D2.dn.x(),         D2.dn.y() - D2.dp, D2.dn.z(),
             0.0,   n.x(),             n.y(),             n.z();

        Vector4d delta = M.inverse() * F;

        double   delta_p = delta[0];
        Vector3d delta_n = {delta[1], delta[2], delta[3]};

        double lambda = 1.0;
        if (n.z() + delta_n.z() < 0.0) {
            // Запрет на выход из диапазона n.z() > 0
            lambda = 0.95 * std::abs(n.z() / delta_n.z());
        }
        delta_p *= lambda;
        delta_n *= lambda;

        n += delta_n; n.normalize();
        p += delta_p;

        double new_ae = get_ae(p, n);

        if (delta_n.norm() < eps_n) {
            break;
        }
        if (std::abs(delta_p) < eps_p) {
            break;
        }
        if (std::abs(ae - new_ae) < eps_ae) {
            ae = new_ae;
            break;
        }
        ae = new_ae;
    }
    return ae;
}

inline double e2f(double ax1, double ax2, double ay1, double ay2) {
    auto [min_x, max_x] = sorted(ax1, ax2);
    auto [min_y, max_y] = sorted(ay1, ay2);
    return 0.5 * (min_x * max_y + max_x * min_y + max_x * max_y - min_x * min_y);
}

using mesh::Side3D;

std::array<double, Side3D::count()> face_fractions(double a_cell, const std::array<double, Side3D::count()>& a_neib) {
    // Ребра вдоль Ox
    double a_bz = edge_fraction(a_cell, a_neib[Side3D::B], a_neib[Side3D::Z]);
    double a_bf = edge_fraction(a_cell, a_neib[Side3D::B], a_neib[Side3D::F]);
    double a_tz = edge_fraction(a_cell, a_neib[Side3D::T], a_neib[Side3D::Z]);
    double a_tf = edge_fraction(a_cell, a_neib[Side3D::T], a_neib[Side3D::F]);

    // Ребра вдоль Oy
    double a_zl = edge_fraction(a_cell, a_neib[Side3D::Z], a_neib[Side3D::L]);
    double a_zr = edge_fraction(a_cell, a_neib[Side3D::Z], a_neib[Side3D::R]);
    double a_fl = edge_fraction(a_cell, a_neib[Side3D::F], a_neib[Side3D::L]);
    double a_fr = edge_fraction(a_cell, a_neib[Side3D::F], a_neib[Side3D::R]);

    // Ребра вдоль Oz
    double a_lb = edge_fraction(a_cell, a_neib[Side3D::L], a_neib[Side3D::B]);
    double a_lt = edge_fraction(a_cell, a_neib[Side3D::L], a_neib[Side3D::T]);
    double a_rb = edge_fraction(a_cell, a_neib[Side3D::R], a_neib[Side3D::B]);
    double a_rt = edge_fraction(a_cell, a_neib[Side3D::R], a_neib[Side3D::T]);

    std::array<double, Side3D::count()> res;
    res[Side3D::L] = e2f(a_lb, a_lt, a_zl, a_fl);
    res[Side3D::R] = e2f(a_rb, a_rt, a_zr, a_fr);
    res[Side3D::B] = e2f(a_lb, a_rb, a_bz, a_bf);
    res[Side3D::T] = e2f(a_lt, a_rt, a_tz, a_tf);
    res[Side3D::Z] = e2f(a_zl, a_zr, a_bz, a_tz);
    res[Side3D::F] = e2f(a_fl, a_fr, a_bf, a_tf);
    return res;
}

} // namespace zephyr::geom