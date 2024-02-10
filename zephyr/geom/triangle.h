#pragma once

#include <array>
#include <vector>

#include <zephyr/geom/vector.h>

namespace zephyr::geom {

/// @brief Представление треугольника. Содержит несколько полезных функций,
/// а также отвечает за отображение из барицентрических координат.
struct Triangle {
protected:
    /// @brief Набор вершин (против часовой для плоскости (x, y))
    std::array<Vector3d, 3> verts;

public:
    /// @brief Конструктор по угловым точкам
    /// @warning Сортировка не осуществляется!
    Triangle(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3);

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..2]
    inline Vector3d &operator[](int idx) {
        return verts[idx];
    }

    /// @brief Прямой доступ к данным по индексу
    /// @param idx in [0..2]
    inline const Vector3d &operator[](int idx) const {
        return verts[idx];
    }

    /// @brief Отображение по барицентрическим координатам
    /// @details Нормировка x1, x2, x3 не требуется
    Vector3d operator()(double x1, double x2, double x3) const;

    /// @brief Отображение по барицентрическим координатам
    /// @details Нормировка x1, x2, x3 не требуется
    Vector3d get(double x1, double x2, double x3) const;

    /// @brief Нормаль к треугольнику (по правилу правой руки)
    Vector3d normal() const;

    /// @brief Нормаль к треугольнику, направлена от точки c
    Vector3d normal(const Vector3d& c) const;

    /// @brief Центр треугольника
    Vector3d center() const;

    /// @brief Площадь треугольника
    double area() const;

    /// @brief Барицентр = центру
    Vector3d centroid() const;

    /// @brief Посчитать объемную долю, которая отсекается от ячейки некоторым
    /// телом, точки которого определяются характеристической функцией inside
    /// @param inside Характеристическая функция: true, если точка находится
    /// внутри тела, иначе false
    /// @param n_points Число тестовых точек, погрешность ~ 1/N.
    double volume_fraction(const std::function<bool(const Vector3d&)>& inside,
                           int n_points = 10000) const;

    /// @brief Интеграл скалярной функции по треугольному элементу
    /// @param n Число подъячеек по осям
    /// @details Сумма по барицентрам 2-го порядка (low accuracy order)
    double integrate_low(const std::function<double(const Vector3d&)>& func, int n) const;

    /// @brief Интеграл скалярной функции по треугольному элементу
    /// @param n Разбиение по сторонам
    /// @details Формула 3-го порядка по 6 узлам (middle accuracy order)
    double integrate_mid(const std::function<double(const Vector3d&)>& func, int n) const;

    /// @brief Интеграл скалярной функции по треугольному элементу
    /// @param n Разбиение по сторонам
    /// @details Формула 6-го порядка по 12 узлам (high accuracy order)
    double integrate_high(const std::function<double(const Vector3d&)>& func, int n) const;

    /// @brief Интеграл скалярной функции по треугольному элементу
    /// @param n Разбиение по сторонам
    /// @details Формула 13-го порядка по 37 узлам (extra-high  accuracy order)
    double integrate_extra(const std::function<double(const Vector3d&)>& func, int n) const;

private:
    /// @brief Отображение по барицентрическим координатам
    /// @details Требуется предворительная нормировка sum x_i = 1
    Vector3d pget(double x1, double x2, double x3) const;
};

} // namespace zephyr::geom