#pragma once

#include <array>

namespace zephyr {
namespace phys {

/// @brief Вектор массовых или объемных концентраций
struct Fractions {
    /// @brief Максимальное число компонент
    static const int max_size = 5;

    /// @brief Конструктор по умолчанию
    Fractions();

    /// @brief Конструктор со списком инициализации
    Fractions(std::initializer_list<double> list);

    /// @brief Конструктор из вектора
    Fractions(const std::vector<double> &list);

    /// @brief Содержит компоненту с индексом idx
    [[nodiscard]] bool has(int idx) const;

    /// @brief Оператор доступа по индексу
    double &operator[](int idx);

    /// @brief Оператор доступа по индексу
    const double &operator[](int idx) const;

    /// @brief Оператор доступа по индексу
    double &operator[](size_t idx);

    /// @brief Оператор доступа по индексу
    const double &operator[](size_t idx) const;

    /// @brief Возвращает true (чистое вещество), если массив содержит
    /// единственную концентрацию больше нуля, в остальных случаях false.
    [[nodiscard]] bool is_pure() const;

    /// @brief Если массив содержит единственную концентрацию,
    /// отличную от нуля, тогда возвращается индекс ненулевого элемента.
    /// Иначе значение -1.
    [[nodiscard]] int index() const;

    /// @brief Нормализовать концентрации (сумма равна единице)
    void normalize();

    /// @brief Обрезать маленькие (< eps) и близкие к единице ( > 1 - eps)
    /// концентрации, затем нормализовать концентрации
    void cutoff(double eps = 1.0e-6);

    [[nodiscard]] size_t get_begin_size() const;

private:
    size_t m_begin_size = max_size;
    std::array<double, max_size> m_data{};
};

std::ostream &operator<<(std::ostream &os, const Fractions &frac);

} // namespace phys
} // namespace zephyr