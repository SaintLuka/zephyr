#pragma once

#include <array>

namespace zephyr { namespace phys {

/// @brief Вектор массовых или объемныых концентраций
struct Fractions {
    /// @brief Максимальное число компонент
    static const int max_size = 5;

    /// @brief Конструктор по умолчанию
    Fractions();

    /// @brief Конструктор со списком инициализации
    Fractions(std::initializer_list<double> list);

    /// @brief Содержит компоненту с индексом idx
    bool has(int idx) const;

    /// @brief Оператор доступа по индексу
    double &operator[](int idx);

    /// @brief Оператор доступа по индексу
    const double &operator[](int idx) const;

    /// @brief Возвращает true (чистое вещество), если массив содержит
    /// единственную концентрацию больше нуля, в остальных случаях false.
    bool is_pure() const;

    /// @brief Если массив содержит единственную концентрацию,
    /// отличную от нуля, тогда возвращается индекс ненулевого элемента.
    /// Иначе значение -1.
    int index() const;

    /// @brief Нормализовать концентрации (сумма равна единице)
    void normalize();

    /// @brief Обрезать маленькие (< eps) и близкие к единице ( > 1 - eps)
    /// концентрации, затем нормализовать концентрации
    void cutoff(double eps = 1.0e-6);

private:
    std::array<double, max_size> m_data;
};

std::ostream& operator<<(std::ostream& os, const Fractions& frac);

} // namespace phys
} // namespace zephyr