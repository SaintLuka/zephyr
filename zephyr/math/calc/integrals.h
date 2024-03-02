#pragma once

#include <functional>

namespace zephyr::math {

enum class IntMethod {
    TRAPZ,
    SIMPS,
    GAUSS
};

// Одномерный интеграл
/// n - число отрезков
double integrate_trapz(const std::function<double(double)>& f,
        double a, double b, int n);

double integrate_simpson(const std::function<double(double)>& f,
        double a, double b, int n);

double integrate_gauss(const std::function<double(double)>& f,
        double a, double b, int n);

double integrate(const std::function<double(double)>& f,
        double a, double b, int n, IntMethod method = IntMethod::SIMPS);

} // namespace zephyr::math