#pragma once

#include <array>
#include <span>
#include <zephyr/geom/vector.h>

namespace zephyr::phys {

class ScalarSet;

/// @brief Вектор массовых или объемных концентраций
class Fractions {
public:
    /// @brief Меньшие объемные/массовые доли не различимы,
    /// при вычитании из единицы остается единица. Около 5.5e-17
    static constexpr double minimal = std::numeric_limits<double>::epsilon() / 4;

    /// @brief Максимальное число компонент
    static constexpr int max_size = 4;

    /// @brief Инициализация нулями
    Fractions() = default;

    /// @brief Конструктор со списком инициализации
    Fractions(std::initializer_list<double> list);

    /// @brief Конструктор из std::array
    Fractions(const std::array<double, max_size> &arr);

    /// @brief Массив нулевых концентраций
    static constexpr Fractions Zero() { return Fractions{}; }

    /// @brief Единственная ненулевая концентрация
    static constexpr Fractions Pure(int idx) {
        Fractions res = Zero(); res.m_data[idx] = 1.0; return res;
    }

    /// @brief Содержит компоненту с индексом idx
    bool has(int idx) const { return m_data[idx] > 0.0; }

    /// @brief Возвращает true (чистое вещество), если массив содержит
    /// единственную концентрацию больше нуля, в остальных случаях false.
    bool is_pure() const;

    /// @brief Установить чистый материал с номером idx
    void set_pure(int idx);

    /// @brief Если массив содержит единственную концентрацию,
    /// отличную от нуля, тогда возвращается индекс ненулевого элемента.
    /// Иначе значение -1.
    int index() const;

    /// @brief Нормализовать концентрации (сумма равна единице)
    void normalize();

    /// @brief Обрезать маленькие (< eps) и близкие к единице ( > 1 - eps)
    /// концентрации, затем нормализовать концентрации
    void cutoff(double eps = 1.0e-12);

    /// @brief Не содержит концентраций?
    bool empty() const;

    /// @brief Число компонент. Индекс последней ненулевой компоненты плюс один.
    /// beta = {1.0, 0.0, 0.0, 2.0, 0.0, 3.0, 0.0, 0.0, 0.0},
    /// тогда beta.count() равно 6!
    int count() const;

    /// @brief Пара индексов материалов для смешанной ячейки
    /// Если ячейка содержит меньше двух материалов, тогда возвращает {-1, -1}.
    /// Если больше двух, тогда возвращает первые два.
    std::tuple<int, int> pair() const;

    // --------------------- Standard container functions ---------------------

    /// @brief Полное число компонент
    constexpr int size() const { return max_size; }

    /// @brief Оператор доступа по индексу
    double &operator[](int idx) { return m_data[idx]; }

    /// @brief Оператор доступа по индексу
    const double &operator[](int idx) const { return m_data[idx]; }

    /// @brief Итератор на начало
    auto begin() { return m_data.begin(); }

    /// @brief Итератор на конец
    auto end() { return m_data.end(); }

    /// @brief Константный итератор на начало
    auto begin() const { return m_data.cbegin(); }

    /// @brief Константный итератор на конец
    auto end() const { return m_data.cend(); }

    /// @brief Преобразование в std::span
    std::span<const double, max_size> span() const { return m_data; }

    /// @brief Преобразование в std::array
    const std::array<double, max_size>& array() const { return m_data; }

private:
    /// @brief Массив данных
    std::array<double, max_size> m_data{};
};

std::ostream &operator<<(std::ostream &os, const Fractions &frac);


/// @brief Набор скалярных величин по числу компонент
class ScalarSet {
    /// @brief Максимальное число компонент
    static constexpr int max_size = Fractions::max_size;

    /// @brief Массив данных
    std::array<double, max_size> m_data{};

public:
    /// @brief Инициализация нулями
    ScalarSet() = default;

    /// @brief Конструктор со списком инициализации
    ScalarSet(std::initializer_list<double> list);

    /// @brief Конструктор из std::array
    ScalarSet(const std::array<double, max_size>& arr);

    /// @brief Конструктор из Fractions
    ScalarSet(const Fractions& frac);

    /// @brief Массив неопределенных значений
    static constexpr ScalarSet NaN() {
        ScalarSet res; res.m_data.fill(NAN); return res;
    }

    /// @brief Все компоненты кроме одной равны нулю.
    static constexpr ScalarSet Pure(int idx, double value) {
        ScalarSet res; res.m_data[idx] = value; return res;
    }

    /// @brief Все компоненты кроме одной равны NAN.
    static constexpr ScalarSet PureNaN(int idx, double value) {
        ScalarSet res = NaN(); res.m_data[idx] = value; return res;
    }

    // --------------------- Standard container functions ---------------------

    /// @brief Полное число компонент
    constexpr int size() const { return max_size; }

    /// @brief Оператор доступа по индексу
    double &operator[](int idx) { return m_data[idx]; }

    /// @brief Оператор доступа по индексу
    double operator[](int idx) const { return m_data[idx]; }

    /// @brief Итератор на начало
    auto begin() { return m_data.begin(); }

    /// @brief Итератор на конец
    auto end() { return m_data.end(); }

    /// @brief Константный итератор на начало
    auto begin() const { return m_data.cbegin(); }

    /// @brief Константный итератор на конец
    auto end() const { return m_data.cend(); }

    /// @brief Преобразование в std::span
    std::span<const double, max_size> span() const { return m_data; }
};

std::ostream &operator<<(std::ostream &os, const ScalarSet &frac);

using geom::Vector3d;

/// @brief Набор векторов по числу компонент
class VectorSet {
public:
    /// @brief Максимальное число компонент
    static constexpr int max_size = Fractions::max_size;

    /// @brief Инициализация нулями
    VectorSet() = default;

    /// @brief Конструктор из вектора
    VectorSet(std::initializer_list<Vector3d> list);

    // --------------------- Standard container functions ---------------------

    /// @brief Полное число компонент
    constexpr int size() const { return max_size; }

    /// @brief Оператор доступа по индексу
    Vector3d& operator[](int idx) { return m_data[idx]; }

    /// @brief Оператор доступа по индексу
    const Vector3d& operator[](int idx) const { return m_data[idx]; }

    /// @brief Итератор на начало
    auto begin() { return m_data.begin(); }

    /// @brief Итератор на конец
    auto end() { return m_data.end(); }

    /// @brief Константный итератор на начало
    auto begin() const { return m_data.cbegin(); }

    /// @brief Константный итератор на конец
    auto end() const { return m_data.cend(); }

private:
    /// @brief Массив данных
    std::array<Vector3d, max_size> m_data{Vector3d::Zero()};
};

std::ostream &operator<<(std::ostream &os, const VectorSet &arr);

} // namespace zephyr