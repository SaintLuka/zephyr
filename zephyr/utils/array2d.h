#pragma once
#include <zephyr/configuration.h>
#include <vector>
#include <array>

namespace zephyr::utils {

/// @brief Двумерный массив с возможностью индексации отрицательными
/// числами как в python.
template <typename T>
class Array2D {
public:
    /// @brief Пустой массив
    Array2D() = default;

    /// @brief Создать массив заданного размера
    explicit Array2D(const std::array<int, 2> sizes, const T& value = {})
        : m_sizes(sizes) {
        z_assert(sizes[0] > 0 && sizes[1] > 0, "Bad sizes constructor");
        m_data.resize(m_sizes[0] * m_sizes[1], value);
    }

    /// @brief Пустой массив?
    bool empty() const { return m_data.empty(); }

    /// @brief Изменить размер массива
    void resize(const std::array<int, 2> sizes, const T& value = {}) {
        z_assert(sizes[0] > 0 && sizes[1] > 0, "Bad sizes resize");
        m_sizes = sizes;
        m_data.clear();
        m_data.resize(m_sizes[0] * m_sizes[1], value);
    }

    /// @brief Размер массива по первой оси
    int size1() const { return m_sizes[0]; }

    /// @brief Размер массива по второй оси
    int size2() const { return m_sizes[1]; }

    /// @brief Размер массива по оси axis
    int size(int axis) const { return m_sizes[axis]; }

    /// @brief Размеры массива
    std::array<int, 2> sizes() const { return m_sizes; }

    /// @brief Получить значение по ссылке
    T& operator()(int i, int j) {
        z_assert(-m_sizes[0] <= i && i < m_sizes[0], "Bad index i");
        z_assert(-m_sizes[1] <= j && j < m_sizes[1], "Bad index j");
        return m_data[reduce(j, m_sizes[1]) * m_sizes[0] + reduce(i, m_sizes[0])];
    }

    /// @brief Получить значение
    const T& operator()(int i, int j) const {
        z_assert(-m_sizes[0] <= i && i < m_sizes[0], "Bad index i");
        z_assert(-m_sizes[1] <= j && j < m_sizes[1], "Bad index j");
        return m_data[reduce(j, m_sizes[1]) * m_sizes[0] + reduce(i, m_sizes[0])];
    }

private:
    /// @brief Привести индекс к корректному диапазону
    static int reduce(int idx, int size) { return idx < 0 ? idx + size : idx; }

    std::array<int, 2> m_sizes{}; ///< Размеры по двум осям
    std::vector<T> m_data{};      ///< Массив с данными
};

} // namespace zephyr::utils