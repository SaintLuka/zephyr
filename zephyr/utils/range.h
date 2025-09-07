#pragma once

#include <algorithm>

/// @brief Основное пространство имен.
namespace zephyr {

/// @brief Простой класс для представления диапазона чисел.
/// @tparam T Целочисленный тип (int, unsigned int, size_t...),
/// по умолчанию используется тип size_t.
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

    static_assert(std::is_integral<T>::value,
                  "Integral type is required for the range");

public:

    /// @brief Примитивный итератор
    class iterator {
    private:
        T m_val; ///< Значение из диапазона

    public:
        /// @brief Категория random_access
        using iterator_category = std::random_access_iterator_tag;
        using difference_type   = std::make_signed_t<T>;
        using value_type = T;
        using pointer    = T*;
        using reference  = T;

        iterator(const T &val) : m_val(val) {}

        operator T() { return m_val; }

        /// @brief Получить значение
        T operator*() const { return m_val; }

        /// @brief Инкремент
        void operator++() { ++m_val; }

        /// @brief Инкремент
        void operator++(int) { ++m_val; }

        /// @brief Декремент
        void operator--() { --m_val; }

        /// @brief Декремент
        void operator--(int) { --m_val; }

        /// @brief Смещение по индексу
        T operator[](difference_type n) const { return m_val + n; }

        /// @brief Инкремент на шаг step
        void operator+=(size_t step) { m_val += step; }

        /// @brief Инкремент на шаг step
        iterator operator+(size_t val) const { return m_val + val; }

        /// @brief Расстояние между итераторами
        difference_type operator-(const iterator &other) const {
            return m_val - other.m_val;
        }

        /// @brief Сравнение
        bool operator!=(const iterator &other) const {
            return m_val != other.m_val;
        }

        /// @brief Сравнение
        bool operator<(const iterator &other) const {
            return m_val < other.m_val;
        }
    };

    /// @brief Конструктор по конечному элементу
    /// @param end Конечный элемент диапазона (не включительно)
    explicit range(const T &end) : m_begin{}, m_end(end) {}

    /// @brief Конструктор по начальному и конечному элементу
    /// @param begin Начальный элемент диапазона
    /// @param end Конечный элемент диапазона (не включительно)
    range(const T &begin, const T &end)
            : m_begin(begin), m_end(end) {}

    /// @brief Начало диапазона
    iterator begin() const { return m_begin; }

    /// @brief Конечный элемент диапазона (не включительно)
    iterator end() const { return m_end; }
};

/// @brief Аналог enumerate из python, можно использовать в циклах со
/// стандартными контейнерами. for (auto [i, val]: some_array) { ... }
/// @tparam Container Стандартный контейнер (vector, array, ...)
template <typename Container>
class enumerate {
public:
    using cont_iter  = typename Container::iterator;
    using value_type = typename Container::value_type;
    using reference  = typename Container::reference;

    enumerate(Container& iterable) : m_container(iterable) { }

    class iterator {
    public:
        iterator(size_t index, cont_iter it) : index(index), it(it) {}

        std::pair<size_t, reference> operator*() const {
            return {index, *it};
        }

        iterator& operator++() {
            ++index;
            ++it;
            return *this;
        }

        bool operator!=(const iterator& other) const {
            return it != other.it;
        }

    private:
        size_t index;
        cont_iter it;
    };

    iterator begin() {
        return iterator(0, std::begin(m_container));
    }

    iterator end() {
        return iterator(0, std::end(m_container));
    }

private:
    Container& m_container;
};

} // namespace zephyr::utils