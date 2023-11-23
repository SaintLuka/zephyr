#pragma once

#include <array>
#include <vector>
#include <zephyr/math/vectorization.h>
#include <zephyr/geom/vector.h>

namespace zephyr::phys {

struct FractionsFlux;

/// @brief Вектор массовых или объемных концентраций
struct Fractions {
    /// @brief Максимальное число компонент
    static const int max_size = 5;

    /// @brief Конструктор по умолчанию
    Fractions();

    /// @brief Конструктор со списком инициализации
    Fractions(std::initializer_list<double> list);

    /// @brief Конструктор из вектора
    explicit Fractions(const std::vector<double> &vec);

    explicit Fractions(const std::array<double, max_size> &arr);

    explicit Fractions(const FractionsFlux &frac_flux);

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

    [[nodiscard]] std::array<double, Fractions::max_size> get_data() const {
        return m_data;
    }

    /// @brief Если массив содержит единственную концентрацию,
    /// отличную от нуля, тогда возвращается индекс ненулевого элемента.
    /// Иначе значение -1.
    [[nodiscard]] int index() const;

    /// @brief Нормализовать концентрации (сумма равна единице)
    void normalize();

    /// @brief Обрезать маленькие (< eps) и близкие к единице ( > 1 - eps)
    /// концентрации, затем нормализовать концентрации
    void cutoff(double eps = 1.0e-6);

    [[nodiscard]] bool empty() const;

    [[nodiscard]] size_t get_size() const;

    std::array<double, max_size> m_data{};

    VECTORIZE(Fractions)
};

std::ostream &operator<<(std::ostream &os, const Fractions &frac);


/// @brief Вектор потока величины нескольких веществ
struct FractionsFlux {
    /// @brief Конструктор по умолчанию
    FractionsFlux();

    /// @brief Конструктор из Fractions
    explicit FractionsFlux(const Fractions &frac);

    /// @brief Конструктор из вектора
    explicit FractionsFlux(const std::vector<double> &vec);

    /// @brief Содержит компоненту с индексом idx
    [[nodiscard]] bool has(int idx) const;

    template<typename T>
    FractionsFlux &operator*=(const T &c) {
        for (double &v: m_data)
            v *= c;

        return *this;
    }

    template<typename T>
    FractionsFlux &operator/=(const T &c) {
        for (double &v: m_data)
            v /= c;

        return *this;
    }

    [[nodiscard]] size_t get_size() const {
        return m_data.size();
    }

    [[nodiscard]] const std::array<double, Fractions::max_size> &get_data() const {
        return m_data;
    }

    std::array<double, Fractions::max_size> m_data{};

    VECTORIZE(FractionsFlux)
};

std::ostream &operator<<(std::ostream &os, const FractionsFlux &frac);

} // namespace zephyr