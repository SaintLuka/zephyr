#pragma once

#include <zephyr/geom/vector.h>

#include <random>

namespace zephyr::geom {

class Random2D;
class QuasiRandom2D;

// TODO: Написать нормально
class Box {
public:

    Vector3d vmin, vmax;

    /// @brief Инициализирует NAN
    Box();

    Box(const Vector3d& _vmin, const Vector3d& _vmax);

    Vector3d center() const;

    Vector3d size() const;
    
    double diameter() const;

    /// @brief Проверяет совпадение z-координаты
    bool is_2D() const;

    /// @brief Проверяет совпадение z-координаты
    bool is_3D() const;

    /// @brief Площадь для двумерного Box
    double area() const;

    /// @brief Объем для трехмерного Box
    double volume() const;

    // Точка внутри?
    bool inside(const Vector3d& p) const;

    /// @brief Сунуть точку внутрь Box
    Vector3d shove_in(const Vector3d& p) const;

    std::vector<double> outline_x() const;

    std::vector<double> outline_y() const;

    void extend(double margin);

    void extend(double margin_x, double margin_y, double margin_z = 0.0);

    Random2D random2D(int seed = 0) const;

    QuasiRandom2D quasiRandom2D() const;
};

std::ostream& operator<<(std::ostream& os, const Box& box);

/// @brief Случайная двумерная последовательность в прямоугольнике
class Random2D {
public:

    Random2D(const Vector3d& vmin, const Vector3d& vmax, int seed = 13);

    Vector3d get();

private:
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> distr_x;
    std::uniform_real_distribution<double> distr_y;
};

/// @brief Квазислучайная двумерная последовательность
class QuasiRandom2D {
public:

    QuasiRandom2D(const Vector3d& vmin, const Vector3d& size);

    Vector3d get();

private:
    Vector3d shift;
    Vector3d step;
    Vector3d vmin;
    Vector3d size;
};

} // namespace zephyr::geom

