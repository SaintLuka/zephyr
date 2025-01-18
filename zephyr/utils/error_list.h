#pragma once

#include <algorithm>
#include <vector>
#include <cmath>

namespace zephyr::utils {

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
        if (x == 0.0 || y == 0.0) {
            return std::max(std::abs(x), std::abs(y));
        }
        return std::abs((x - y) / (std::abs(x) + std::abs(y)));
    }

    /// @brief Добавить пару значений для сравнения
    void add(const std::pair<double, double>& p) {
        val1.push_back(p.first);
        val2.push_back(p.second);
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

    int size() const {
        return errs.size();
    }

    double operator[](int idx) const {
        return errs[idx];
    }

    static void print_arr(const std::vector<double>& arr) {
        std::cout << "[ ";
        for (int i = 0; i < arr.size() - 1; ++i) {
            std::cout << arr[i] << ", ";
        }
        std::cout << arr.back() << " ]";
    }

    void print(const std::string& tab = "\t") const {
        std::cout << tab << "vals 1: ";
        print_arr(val1);
        std::cout << "\n" << tab << "vals 2: ";
        print_arr(val2);
        std::cout << "\n" << tab << "errors: ";
        print_arr(errs);
        std::cout << "\n";
    }

private:
    ///@brief Массив ошибок
    std::vector<double> val1;
    std::vector<double> val2;
    std::vector<double> errs;
};

}