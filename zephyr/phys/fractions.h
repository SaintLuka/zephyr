#pragma once

#include <array>
#include <vector>
#include <zephyr/math/vectorization.h>
#include <zephyr/geom/vector.h>

namespace zephyr::phys {

struct ScalarSet;

/// @brief Вектор массовых или объемных концентраций
struct Fractions {
    /// @brief Меньшие объемные/массовые доли не различимы,
    /// при вычитании из единицы остается единица. Около 5.5e-17
    constexpr static const double minimal = std::numeric_limits<double>::epsilon() / 4;

    /// @brief Максимальное число компонент
    static const int max_size = 3;

    /// @brief Массив данных
    std::array<double, max_size> m_data{};


    /// @brief Инициализация нулями
    Fractions();

    /// @brief Конструктор со списком инициализации
    Fractions(std::initializer_list<double> list);

    /// @brief Конструктор из std::vector
    explicit Fractions(const std::vector<double> &vec);

    /// @brief Конструктор из std::array
    explicit Fractions(const std::array<double, max_size> &arr);

    /// @brief Конструктор из набора скаляров
    explicit Fractions(const ScalarSet &scalars);

    /// @brief Содержит компоненту с индексом idx
    bool has(int idx) const;

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
    bool is_pure() const;

    /// @brief Установить чистый материал с номером idx
    void set_pure(int idx);

    /// @brief Если массив содержит единственную концентрацию,
    /// отличную от нуля, тогда возвращается индекс ненулевого элемента.
    /// Иначе значение -1.
    int index() const;

    /// @brief Количество компонент.
    /// Индекс последней ненулевой компоненты плюс один.
    /// beta = {1.0, 0.0, 0.0, 2.0, 0.0, 3.0, 0.0, 0.0, 0.0},
    /// тогда beta.count() равно 6!
    int count() const;

    /// @brief Количество ненулевых компонент.
    /// Индекс последней ненулевой компоненты плюс один.
    /// beta = {1.0, 0.0, 0.0, 2.0, 0.0, 3.0, 0.0, 0.0, 0.0},
    /// тогда beta.nonzero() равно 3!
    int nonzero() const;

    /// @brief Пара индексов материалов для смешаной ячейки
    /// Если ячейка содержит меньше двух материалов, тогда возвращает {-1, -1}.
    /// Если больше двух, тогда возвращает первые два.
    std::array<int, 2> pair() const;

    /// @brief Нормализовать концентрации (сумма равна единице)
    void normalize();

    /// @brief Обрезать маленькие (< eps) и близкие к единице ( > 1 - eps)
    /// концентрации, затем нормализовать концентрации
    void cutoff(double eps = 1.0e-12);

    /// @brief Не содержит концентраций
    bool empty() const;

    /// @brief Указатель на данные
    double* data() { return m_data.data(); }

    /// @brief Константный указатель на данные
    const double* data() const { return m_data.data(); }

    /// @brief Ссылка на массив данных
    const std::array<double, Fractions::max_size>& data_ref() const {
        return m_data;
    }

    VECTORIZE(Fractions)
};

std::ostream &operator<<(std::ostream &os, const Fractions &frac);


/// @brief Вектор потока величины нескольких веществ
struct ScalarSet {
    /// @brief Массив данных
    std::array<double, Fractions::max_size> m_data{};


    /// @brief Конструктор по умолчанию
    /// Инициализирует нулями.
    ScalarSet();

    /// @brief Установить значение
    ScalarSet(double val);

    /// @brief Конструктор со списком инициализации
    ScalarSet(std::initializer_list<double> list);

    /// @brief Установить единственное значение по индексу idx,
    /// остальные компоненты инициализируются нулями.
    ScalarSet(double val, int idx);

    /// @brief Конструктор из Fractions
    explicit ScalarSet(const Fractions &frac);

    /// @brief Конструктор из вектора
    explicit ScalarSet(const std::vector<double> &vec);

    /// @brief Оператор доступа по индексу
    double &operator[](size_t idx);

    /// @brief Оператор доступа по индексу
    const double &operator[](size_t idx) const;

    template<typename T>
    ScalarSet &operator*=(const T &c) {
        for (double &v: m_data)
            v *= c;

        return *this;
    }

    template<typename T>
    ScalarSet &operator/=(const T &c) {
        for (double &v: m_data)
            v /= c;

        return *this;
    }

    ScalarSet &operator+=(const ScalarSet &c) {
        for (int i = 0; i < size(); ++i) {
            m_data[i] += c[i];
        }

        return *this;
    }

    /// @brief Ссылка на массив данных
    const std::array<double, Fractions::max_size> &data_ref() const {
        return m_data;
    }

    VECTORIZE(ScalarSet)
};

std::ostream &operator<<(std::ostream &os, const ScalarSet &frac);

} // namespace zephyr