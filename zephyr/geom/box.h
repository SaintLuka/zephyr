#pragma once

#include <zephyr/geom/vector.h>

namespace zephyr::geom {

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

    /// @brief Площадь для двумерного Box
    double area() const;

    /// @brief Объем для трехмерного Box
    double volume() const;

    bool inside(const Vector3d& p) const;

    void extend(double margin);

    void extend(double margin_x, double margin_y, double margin_z = 0.0);

    QuasiRandom2D quasiRandom2D() const;
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

