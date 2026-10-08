#pragma once

#include <memory>
#include <random>

#include <zephyr/geom/box.h>

namespace zephyr::math {

class Random {
public:
    using Ptr = std::shared_ptr<Random>;

    virtual ~Random() = default;

    virtual geom::Vector3d next() = 0;

    /// @brief Равномерное случайное распределение в прямоугольнике/кубоиде
    static Random::Ptr Uniform(const geom::Box& box, int seed = 13);

    /// @brief Квазислучайное распределение в прямоугольнике/кубоиде
    static Random::Ptr Quasi(const geom::Box& box);
};

/// @brief Равномерно распределенные точки в прямоугольнике
class Random2D : public Random {
public:
    Random2D(const geom::Box& box, int seed = 13);

    Random2D(const geom::Vector3d& vmin,
             const geom::Vector3d& vmax, int seed = 13);

    ~Random2D() override = default;

    geom::Vector3d next() final;

private:
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> distr_x;
    std::uniform_real_distribution<double> distr_y;
};

/// @brief Равномерно распределенные точки в кубоиде
class Random3D : public Random {
public:
    Random3D(const geom::Box& box, int seed = 13);

    Random3D(const geom::Vector3d& vmin,
             const geom::Vector3d& vmax, int seed = 13);

    ~Random3D() override = default;

    geom::Vector3d next() final;

private:
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> distr_x;
    std::uniform_real_distribution<double> distr_y;
    std::uniform_real_distribution<double> distr_z;
};

/// @brief Квазислучайная последовательность точек в прямоугольнике
class QuasiRandom2D : public Random {
public:
    QuasiRandom2D(const geom::Box& box);

    QuasiRandom2D(const geom::Vector3d& vmin,
                  const geom::Vector3d& size);

    ~QuasiRandom2D() override = default;

    /// @brief Получить следующую точку
    geom::Vector3d next() final;

private:
    geom::Vector3d shift;
    geom::Vector3d step;
    geom::Vector3d vmin;
    geom::Vector3d size;
};

} // namespace zephyr::math