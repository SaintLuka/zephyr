#pragma once

namespace zephyr::utils {

/// @class Простой класс для представления диапазона чисел.
/// По умолчанию используется для целочисленного типа int, можно использовать
/// для генерации чисел из size_t. Вообще говоря, класс нужен, чтобы кидать
/// на вход функции threads::for_each().
/// Пример, вывести числа от 0 до 9 (включительно):
/// @code
/// for (auto i: range(10)) {
///     std::cout << i << "\n";
/// }
/// @endcode
template<typename T = size_t>
class range {
private:
    T m_begin;  ///< Начальный элемент диапазона
    T m_end;    ///< Конечный элемент диапазона (не включительно)

    /// @brief В качестве шаблонного параметра класса должен использоваться
    /// целочисленный тип (int, unsigned int, size_t...)
    static_assert(std::is_integral<T>::value,
                  "Integral type is required for the range");

public:

    /// @brief Примитивный итератор
    class iterator {
    private:
        T m_val; ///< Значение из диапазона

    public:
        iterator(const T &val) : m_val(val) {}

        inline operator T() { return m_val; }

        /// @brief Получить значение
        const T &operator*() const { return m_val; }

        /// @brief Инкремент
        inline void operator++() { ++m_val; }

        /// @brief Инкремент
        inline void operator++(int) { ++m_val; }

        /// @brief Инкремент на шаг step
        inline void operator+=(size_t step) { m_val += step; }

        /// @brief Инкремент на шаг step
        inline iterator operator+(size_t val) const {
            return m_val + val;
        }

        /// @brief Расстояние между итераторами
        inline int operator-(const iterator &other) const {
            return int(m_val - other.m_val);
        }

        /// @brief Сравнение
        inline bool operator!=(const iterator &other) {
            return m_val != other.m_val;
        }

        /// @brief Сравнение
        inline bool operator<(const iterator &other) {
            return m_val < other.m_val;
        }
    };

    /// @brief Конструктор по конечному элементу
    /// @param end Конечный элемент диапазона (не включительно)
    explicit range(const T &end)
            : m_begin{}, m_end(end) {}

    /// @brief Конструктор по начальному и конечному элементу
    /// @param begin Начальный элемент диапазона
    /// @param end Конечный элемент диапазона (не включительно)
    range(const T &begin, const T &end)
            : m_begin(begin), m_end(end) {}

    /// @brief Начало диапазона
    inline iterator begin() const {
        return m_begin;
    }

    /// @brief Конечный элемент диапазона (не включительно)
    inline iterator end() const {
        return m_end;
    }
};

} // namespace zephyr::utils