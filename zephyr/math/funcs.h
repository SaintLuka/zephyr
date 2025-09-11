/// @file Некоторые математические функции
/// Параметрические функции часто в виде функторов

#include <cmath>
#include <tuple>
#include <vector>

namespace zephyr::math {

/// @brief Минимум из трёх величин
inline double min(double x, double y, double z) {
    return std::min(x, std::min(y, z));
}

/// @brief Максимум из трёх величин
inline double max(double x, double y, double z) {
    return std::max(x, std::max(y, z));
}

/// @brief Функция знака
inline double sign(double x) {
    return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0);
}

/// @brief Функция Хевисайда
/// @details Доопределена нулем в нуле.
inline double heav(double x) {
    return x > 0.0 ? 1.0 : 0.0;
}

/// @brief Возвращает упорядоченную пару (минимум/максимум)
inline std::tuple<double, double> sorted(double a, double b) {
    if (a < b) return {a, b};
    else       return {b, a};
}

/// @brief Возвращает упорядоченную тройку
inline std::tuple<double, double, double> sorted(double a, double b, double c) {
    if (a > b) std::swap(a, b);
    if (a > c) std::swap(a, c);
    if (b > c) std::swap(b, c);
    return {a, b, c};
}

/// @brief Ограничить значение x сверху и снизу.
/// Возвращает число между min(x1, x2) и max(x1, x2).
inline double between(double x, double x1, double x2) {
    auto[x_min, x_max] = sorted(x1, x2);
    return std::max(x_min, std::min(x, x_max));
}

/// @brief Гладкая функция Хевисайда с синусом
/// @param w Ширина (>= 0.0)
inline double heav_s(double x, double w = 1.0) {
    if (x <= 0.0) return 0.0;
    if (x >= w)   return 1.0;
    return std::sin(M_PI_2 * x / w);
}

/// @brief Гладкая функция Хевисайда с полиномом
/// @param w Ширина (>= 0.0)
template <int n = 4>
inline double heav_p(double x, double w = 1.0) {
    if (x <= 0.0) return 0.0;
    if (x >= w)   return 1.0;
    return 1.0 - std::pow(1.0 - x / w, n);
}

/// @brief Гладкая функция Хевисайда с корнем
/// @param w Ширина (>= 0.0)
template <int n = 3>
inline double heav_r(double x, double w = 1.0) {
    if (x <= 0.0) return 0.0;
    if (x >= w)   return 1.0;
    double y = n == 3 ? std::cbrt(x / w) :
              (n == 2 ? std::sqrt(x / w) : std::pow(x, 1.0 / n));
    return y * (n + 1 - x / w) / n;
}

/// @brief Гладкая сигмоида с синусом
/// @param w Ширина (>= 0.0)
inline double sign_s(double x, double w = 1.0) {
    if (std::abs(x) >= w) return sign(x);
    return std::sin(M_PI_2 * x / w);
}

/// @brief Гладкая функция Хевисайда с полиномом
/// @param w Ширина (>= 0.0)
template <int n = 4>
inline double sign_p(double x, double w = 1.0) {
    if (std::abs(x) >= w) return sign(x);
    return sign(x) * (1.0 - std::pow(1.0 - std::abs(x / w), n));
}

/// @brief Гладкая функция Хевисайда с корнем
/// @param w Ширина (>= 0.0)
template <int n = 3>
inline double sign_r(double x, double w = 1.0) {
    if (std::abs(x) >= w) return sign(x);

    double xi = std::abs(x / w);
    double y = n == 3 ? std::cbrt(xi) :
              (n == 2 ? std::sqrt(xi) : std::pow(xi, 1.0 / n));
    return sign(x) * y * (n + 1 - xi) / n;
}

} // namespace zephyr::math