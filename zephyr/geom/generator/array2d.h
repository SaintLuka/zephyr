#pragma once

#include <array>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <zephyr/configuration.h>

namespace zephyr::geom::generator {

/// @brief Две оси четырёхугольного блока
enum class Axis : int { X = 0, Y = 1 };

/// @brief Стороны четырёхугольного блока
enum class Side : int {
    L = 0, LEFT   = 0,
    R = 1, RIGHT  = 1,
    B = 2, BOTTOM = 2,
    T = 3, TOP    = 3,
};

/// @brief Список всех сторон четырёхугольника
constexpr std::array sides_2D = {Side::L, Side::R, Side::B, Side::T};

// Выбрать ось блока, вдоль которой лежит сторона
constexpr Axis to_axis(Side side) {
    return (side == Side::B || side == Side::T) ? Axis::X : Axis::Y;
}

/// @brief Пара величин, присвоенная осям блока
template <typename T>
class AxisPair {
public:
    /// @brief Конструктор по умолчанию
    AxisPair() = default;

    /// @brief Конструктор от initializer_list
    constexpr AxisPair(std::initializer_list<T> list) {
        if (list.size() != 2)
            throw std::invalid_argument("Pair requires exactly 2 elements");
        std::copy_n(list.begin(), 2, values.begin());
    }

    /// @brief Конструктор от std::array
    constexpr AxisPair(const std::array<T, 2>& arr) : values(arr) {}

    /// @brief Неявное приведение к std::array
    operator const std::array<T, 2>&() const { return values; }

    /// @brief Доступ по оси
    T& operator[](Axis axis) {
        return values[static_cast<int>(axis)];
    }

    /// @brief Доступ по оси
    const T& operator[](Axis axis) const {
        return values[static_cast<int>(axis)];
    }

    /// @brief Доступ по стороне
    T& operator[](Side side) {
        return values[static_cast<int>(to_axis(side))];
    }

    /// @brief Доступ по стороне
    const T& operator[](Side side) const {
        return values[static_cast<int>(to_axis(side))];
    }

private:
    /// @brief Пара значений
    std::array<T, 2> values{};
};

/// @brief Двумерный массив с возможностью индексации отрицательными
/// индексами как в python. Используется для таблиц узлов.
template <typename T>
class Array2D {
public:
    /// @brief Пустой массив
    Array2D() = default;

    /// @brief Создать массив заданного размера
    explicit Array2D(const std::array<int, 2> sizes, const T& value)
        : m_sizes(sizes) {
        z_assert(sizes[0] > 0 && sizes[1] > 0, "Bad sizes constructor");
        m_data.resize(m_sizes[0] * m_sizes[1], value);
    }

    /// @brief Создать массив из существующего буфера
    explicit Array2D(const std::array<int, 2> sizes, std::vector<T>&& data)
        : m_sizes(sizes) {
        z_assert(sizes[0] > 0 && sizes[1] > 0, "Bad sizes constructor");
        if (sizes[0] * sizes[1] != data.size()) {
            throw std::runtime_error("Array2D: bad sizes for buffer");
        }
        m_data = std::move(data);
    }

    /// @brief Пустой массив?
    bool empty() const { return m_data.empty(); }

    /// @brief Очистить массив
    void clear() {
        m_sizes = {0, 0};
        m_data.clear();
    }

    /// @brief Изменение размеров массива (выполняется с полной очисткой)
    void resize(const std::array<int, 2> sizes) {
        z_assert(sizes[0] > 0 && sizes[1] > 0, "Bad sizes resize");
        m_sizes = sizes;
        m_data.clear();
        m_data.resize(m_sizes[0] * m_sizes[1]);
    }

    /// @brief Изменение размеров массива (выполняется с полной очисткой)
    void resize(const std::array<int, 2> sizes, const T& value) {
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
    int size(Axis axis) const { return m_sizes[static_cast<int>(axis)]; }

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

    /// @brief Угловая вершина блока
    T& corner(int v_idx) {
        switch (v_idx) {
            case 0: return operator()( 0,  0);
            case 1: return operator()(-1,  0);
            case 2: return operator()( 0, -1);
            case 3: return operator()(-1, -1);
            default: throw std::runtime_error("Invalid base node index");
        }
    }

    /// @brief Угловая вершина блока
    const T& corner(int v_idx) const {
        switch (v_idx) {
            case 0: return operator()( 0,  0);
            case 1: return operator()(-1,  0);
            case 2: return operator()( 0, -1);
            case 3: return operator()(-1, -1);
            default: throw std::runtime_error("Invalid base node index");
        }
    }

    /// @brief Граничная вершина блока
    T& boundary(Side side, int idx) {
        switch (side) {
            case Side::B: return operator()(idx, 0);
            case Side::R: return operator()(-1, idx);
            case Side::T: return operator()(-idx - 1, -1);
            case Side::L: return operator()(0, -idx - 1);
            default: throw std::runtime_error("Invalid base node index");
        }
    }

    /// @brief Вершина со второго ряда от границы
    T& near_boundary(Side side, int idx) {
        switch (side) {
            case Side::B: return operator()(idx, 1);
            case Side::R: return operator()(-2, idx);
            case Side::T: return operator()(-idx - 1, -2);
            case Side::L: return operator()(1, -idx - 1);
            default: throw std::runtime_error("Invalid base node index");
        }
    }

    /// @brief Указатель на данные
    const T* data() const { return m_data.data(); }

private:
    /// @brief Привести индекс к корректному диапазону
    static int reduce(int idx, int size) { return idx < 0 ? idx + size : idx; }

    std::array<int, 2> m_sizes{}; ///< Размеры по двум осям
    std::vector<T>     m_data{};  ///< Массив с данными
};

} // zephyr::geom::generator