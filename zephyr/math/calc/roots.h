#pragma once

#include <cmath>
#include <functional>
#include <stdexcept>
#include <optional>
#include <iomanip>
#include <complex>
#include <array>

namespace zephyr::math {

/// @brief Найти корень уравнения func(x) = 0.0
/// @param tol Относительная погрешность
inline double dichotomy(const std::function<double(double)>& func,
                        double x_min, double x_max, double tol = 1.0e-3) {
    double eps = tol * (x_max - x_min);

    double f_min = func(x_min);
    double f_max = func(x_max);

    double delta = tol * std::abs(f_max - f_min);
    if (std::abs(f_min) < delta) { return x_min; }
    if (std::abs(f_max) < delta) { return x_max; }

    if (f_min * f_max > 0.0) {
        std::cerr << "dichotomy: f_min = " << f_min << " " << f_max << std::endl;
        throw std::runtime_error("dichotomy: wrong assumption");
    }

    while (x_max - x_min > eps) {
        double x = 0.5 * (x_min + x_max);
        double f_avg = func(x);

        if (f_min * f_avg <= 0.0) {
            x_max = x;
            f_max = f_avg;
        } else {
            x_min = x;
            f_min = f_avg;
        }
    }

    return 0.5 * (x_min + x_max);
}

struct NewtonOpts {
    double tol = 1e-8;   // Погрешность по невязке
    double atol = 1e-8;  // Погрешность по аргументу
    int max_iters = 30;  // Максимальное количество итераций

    std::optional<double> x_min = std::nullopt;  // Минимальное значение аргумента
    std::optional<double> x_max = std::nullopt;  // Максимальное значение аргумента
};

template<bool verbose = false, typename Func, typename Deriv>
double newton_method(Func f, Deriv df, double x0, const NewtonOpts& options = {}) {
    double x = x0;
    double x_prev;
    double fx = f(x);
    double dfx;

    // Проверка начального приближения
    if (options.x_min && x < options.x_min.value()) {
        throw std::invalid_argument("Начальное приближение меньше минимальной границы");
    }
    if (options.x_max && x > options.x_max.value()) {
        throw std::invalid_argument("Начальное приближение больше максимальной границы");
    }

    if constexpr (verbose) {
        std::cout << "Метод Ньютона:\n";
        std::cout << std::setw(5) << "Iter" << std::setw(15) << "x"
                  << std::setw(15) << "f(x)" << std::setw(15) << "f'(x)"
                  << std::setw(15) << "Δx" << "\n";
        std::cout << std::string(65, '-') << "\n";
    }

    size_t iter = 0;
    bool converged = false;

    for (; iter < options.max_iters; ++iter) {
        x_prev = x;
        dfx = df(x);

        // Итерация Ньютона
        x -= fx / dfx;

        // Проверка границ
        if (options.x_min && x < options.x_min.value()) {
            x = options.x_min.value();
        }
        if (options.x_max && x > options.x_max.value()) {
            x = options.x_max.value();
        }

        fx = f(x);
        double delta_x = std::abs(x - x_prev);

        if constexpr (verbose) {
            std::cout << std::setw(5) << iter << std::setw(15) << x
                      << std::setw(15) << fx << std::setw(15) << dfx
                      << std::setw(15) << delta_x << "\n";
        }

        // Критерии остановки
        if (std::abs(fx) < options.tol || delta_x < options.tol) {
            converged = true;
            break;
        }
    }

    if constexpr (verbose) {
        if (iter == options.max_iters && !converged) {
            std::cerr << "Превышено максимальное количество итераций\n";
        }
    }

    return x;
}

using complex = std::complex<double>;
using roots_t = std::array<complex, 3>;

static constexpr double eps = 1e-14;

// Проверка на ноль с учетом погрешности
inline bool is_zero(double x) {
    return std::abs(x) < eps;
}

// Проверка на вещественное число
inline bool is_real(const complex& z) {
    return std::abs(z.imag()) < eps;
}

// Нормализация вещественного числа (убираем мнимую часть если она близка к нулю)
inline complex normalized(const complex& z) {
    return is_real(z) ? complex{z.real(), 0.0} : z;
}

// Решение канонического кубического уравнения: y³ + py + q = 0
inline roots_t cardano(double p, double q) {
    // Дискриминант
    double Q = q * q / 4.0 + p * p * p / 27.0;

    if (is_zero(Q)) {
        // Кратные корни
        if (is_zero(p) && is_zero(q)) {
            // Тройной корень
            // std::cout << "case 1\n";
            return {complex{}, complex{}, complex{}};
        } else {
            // Двойной корень и простой
            // std::cout << "case 2\n";
            double u = std::cbrt(-q / 2.0);
            complex double_root = {2.0 * u, 0.0};
            complex single_root = {-u, 0.0};
            return {double_root, double_root, single_root};
        }
    } else if (Q > 0) {
        // Один вещественный и два комплексно-сопряженных корня
        // std::cout << "case 3\n";
        double sqrt_Q = std::sqrt(Q);
        double u = std::cbrt(-q / 2.0 + sqrt_Q);
        double v = std::cbrt(-q / 2.0 - sqrt_Q);

        complex real_root = complex{u + v, 0.0};
        complex complex_root1 = complex{-(u + v) / 2.0,  (u - v) * std::sqrt(3.0) / 2.0};
        complex complex_root2 = complex{-(u + v) / 2.0, -(u - v) * std::sqrt(3.0) / 2.0};

        return {real_root, normalized(complex_root1), normalized(complex_root2)};
    } else {
        // Три вещественных корня (неприводимый случай)
        // std::cout << "case 4\n";
        double r = std::sqrt(-p * p * p / 27.0);
        double phi = std::acos(-q / (2.0 * r));

        complex root1 = {2.0 * std::cbrt(r) * std::cos(phi / 3.0), 0.0};
        complex root2 = {2.0 * std::cbrt(r) * std::cos((phi + 2.0 * M_PI) / 3.0), 0.0};
        complex root3 = {2.0 * std::cbrt(r) * std::cos((phi + 4.0 * M_PI) / 3.0), 0.0};

        return {root1, root2, root3};
    }
}

// Решение кубического уравнения: ax³ + bx² + cx + d = 0
inline roots_t solve_cubic(double a, double b, double c, double d) {
    // Нормализуем коэффициенты
    if (is_zero(a)) {
        throw std::invalid_argument("Коэффициент a не может быть нулевым");
    }

    b /= a;
    c /= a;
    d /= a;

    // Приводим к виду: y³ + py + q = 0. Подстановка: x = y - b/3
    double p = c - b * b / 3.0;
    double q = 2.0 * b * b * b / 27.0 - b * c / 3.0 + d;

    auto roots = cardano(p, q);
    roots[0] -= b / 3.0;
    roots[1] -= b / 3.0;
    roots[2] -= b / 3.0;

    return roots;
}

// Проверка решения подстановкой
inline void verify_cubic_roots(double a, double b, double c, double d, const roots_t& roots) {
    std::cout << "\nПроверка решения:\n";
    for (const auto& root : roots) {
        complex value = a * root * root * root +
                      b * root * root +
                      c * root +
                      d;
        std::cout << "f(" << root << ") = " << value;
        if (std::abs(value) < eps) {
            std::cout << " ✓" << std::endl;
        } else {
            std::cout << " ✗" << std::endl;
        }
    }
}


} // namespace zephyr::math