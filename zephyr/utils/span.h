#pragma once

#include <iostream>
#include <vector>
#include <array>

namespace zephyr::utils {

/// @brief View для последовательно расположенных данных
template<typename T>
class span {
public:
    /// @brief Конструктор
    span(T* ptr, size_t size) : m_ptr(ptr), m_size(size) { }

    /// @brief Доступ к элементу
    T& operator[](size_t index) const { return m_ptr[index]; }

    /// @brief Длина диапазона
    size_t size() const { return m_size; }

    /// @brief Неявное преобразование в std::vector<T>
    operator std::vector<std::remove_const_t<T>>() const {
        std::vector<std::remove_const_t<T>> result(m_size);
        for (size_t i = 0; i < m_size; ++i) { result[i] = m_ptr[i]; }
        return result;
    }

    /// @brief Неявное преобразование в std::array<T, N>
    template<size_t N>
    operator std::array<std::remove_const_t<T>, N>() const {
        std::array<std::remove_const_t<T>, N> result{};
        const size_t max_i = std::min(m_size, N);
        for (size_t i = 0; i < max_i; ++i) { result[i] = m_ptr[i]; }
        return result;
    }

    /// @brief Операция присваивания из std::array
    template<size_t N>
    void operator=(const std::array<T, N>& arr) const {
        const size_t max_i = std::min(m_size, N);
        for (size_t i = 0; i < max_i; ++i) { m_ptr[i] = arr[i]; }
    }

    /// @brief Операция присваивания из std::vector
    void operator=(const std::vector<T>& arr) const {
        const size_t max_i = std::min(m_size, size_t{arr.size()});
        for (size_t i = 0; i < max_i; ++i) { m_ptr[i] = arr[i]; }
    }

private:
    T*     m_ptr;  ///< Указатель на начало данных
    size_t m_size; ///< Число элементов
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const span<T>& sp) {
    os << "{ ";
    for (size_t i = 0; i < sp.size() - 1; ++i) {
        os << sp[i] << ", ";
    }
    os << sp[sp.size() - 1] << " }";
    return os;
}

} // namespace zephyr::utils