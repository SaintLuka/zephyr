#pragma once

#include <algorithm>
#include <vector>
#include <cmath>

namespace zephyr { namespace utils {

/// @brief Небольшой класс для проведения тестов.
/// Содержит массив относительных ошибок при сравнении пар значений.
struct ErrorList {
public:
    /// @brief Конструктор по умолчанию
    ErrorList() = default;

    /// @brief Конструктор по списку
    ErrorList(const std::initializer_list<std::pair<double, double>>& list) {
        for (auto& p: list) {
            add(p);
        }
    }

    /// @brief Относительная погрешность при сравнении двух чисел
    static double error(double x, double y) {
        if (x == 0.0 && y == 0.0) {
            return 0.0;
        }
        return std::abs((x - y) / (std::abs(x) + std::abs(y)));
    }

    /// @brief Добавить пару значений для сравнения
    void add(const std::pair<double, double>& p) {
        errs.push_back(error(p.first, p.second));
    }

    /// @brief Добавить пару значений для сравнения
    void operator +=(const std::pair<double, double>& p) {
        add(p);
    }

    /// @brief Максимальная погрешность
    double max() const {
        return *std::max_element(errs.begin(), errs.end());
    }

    /// @brief Погрешности не удовлетворяют условиям
    bool is_ok(double tolerance) {
        // Содержит NAN
        for (auto& v: errs) {
            if (std::isnan(v)) {
                return false;
            }
        }

        // Максимальная ошибка больше tolerance
        return max() < tolerance;
    }

private:
    ///@brief Массив ошибок
    std::vector<double> errs;
};

}
}