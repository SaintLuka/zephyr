#pragma once

#include <zephyr/geom/vector.h>

namespace zephyr::math {
class Random2D;
class QuasiRandom2D;
}

namespace zephyr::geom {

/// @brief Ограничивающий кубоид (bounding box).
/// @details Очень полезная структура для хранения ограничивающего
/// прямоугольника или кубоида (bounding box). Содержит пару полезных функций.
/// Ящик считается двумерным прямоугольником, если vmin.z() == vmax.z(),
/// условие справедливо, в том числе для бесконечных и пустых ящиков.
/// В обратном случае vmin.z() != vmax.z() ящик считается трёхмерном.
struct Box final {
    Vector3d vmin = {NAN, NAN, NAN};
    Vector3d vmax = {NAN, NAN, NAN};

    /// @brief Статический конструктор с неопределенными значениями
    static Box NaN();

    /// @brief Статический конструктор с нулевыми значениями
    static Box Zero();

    /// @brief Минимально возможный bounding box, предполагает
    /// последовательное увеличение путем помещения новых точек.
    /// vmin = {+inf, +inf, +inf}
    /// vmax = {-inf, -inf, -inf}
    /// @param dim Размерность ящика (2 или 3)
    static Box Empty(int dim);

    /// @brief Bounding box охватывающий всё пространство
    /// vmin = {-inf, -inf, -inf}
    /// vmax = {+inf, +inf, +inf}
    /// @param dim Размерность ящика (2 или 3)
    static Box Infinite(int dim);

    /// @brief Центр ящика
    Vector3d center() const;

    /// @brief Вектор размеров ящика по осям координат
    Vector3d size() const;

    /// @brief Расстояние между крайними точками ящика
    double diameter() const;

    /// @brief Проверяет совпадение z-координаты
    bool is_2D() const;

    /// @brief Проверяет совпадение z-координаты
    bool is_3D() const;

    /// @brief Площадь для двумерного Box
    double area() const;

    /// @brief Объем для трехмерного Box
    double volume() const;

    /// @brief Точка находится внутри ящика?
    bool inside(const Vector3d& p) const;

    /// @brief Ограничить координаты точки до границ ящика
    Vector3d shove_in(const Vector3d& p) const;

    /// @brief Расширить ящик, чтобы захватить точку
    void capture(const Vector3d& p);

    /// @brief Замкнутая линия x-координат границ бокса
    std::vector<double> outline_x() const;

    /// @brief Замкнутая линия y-координат границ бокса
    std::vector<double> outline_y() const;

    /// @brief Расширить границы на долю margin по каждой координате
    void extend(double margin);

    /// @brief Расширить границы на долю margin по каждой координате
    void extend(double margin_x, double margin_y, double margin_z = 0.0);

    /// @brief Создать генератор случайных чисел внутри прямоугольника
    math::Random2D random2D(int seed = 0) const;

    /// @brief Создать генератор квазислучайной последовательности внутри
    /// прямоугольника
    math::QuasiRandom2D quasiRandom2D() const;
};

/// @brief Вывод ящика в консоль
std::ostream& operator<<(std::ostream& os, const Box& box);

} // namespace zephyr::geom

