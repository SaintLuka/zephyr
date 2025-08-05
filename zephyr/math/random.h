#pragma once

#include <random>

#include <zephyr/geom/vector.h>

namespace zephyr::geom {}

namespace zephyr::math {

/// @brief Равномерно распределенные точки в прямоугольнике
class Random2D {
public:
    Random2D(const geom::Vector3d& vmin,
             const geom::Vector3d& vmax, int seed = 13);

    geom::Vector3d get();

private:
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> distr_x;
    std::uniform_real_distribution<double> distr_y;
};

/// @brief Квазислучайная последовательность точек в прямоугольнике
class QuasiRandom2D {
public:

    QuasiRandom2D(const geom::Vector3d& vmin,
                  const geom::Vector3d& size);

    /// @brief Получить следующую точку
    geom::Vector3d get();

private:
    geom::Vector3d shift;
    geom::Vector3d step;
    geom::Vector3d vmin;
    geom::Vector3d size;
};

} // namespace zephyr::math