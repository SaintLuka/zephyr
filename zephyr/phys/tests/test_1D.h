#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/matter/eos/eos.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;

/// @class Абстрактный класс одномерного теста
class Test1D {
public:
    /// @brief Получить название теста
    virtual std::string name() const = 0;

    /// @brief Левая граница области
    virtual double xmin() const = 0;

    /// @brief Правая граница области
    virtual double xmax() const = 0;

    /// @brief Конечный момент времени
    virtual double max_time() const = 0;

    /// @brief Получить положение разрыва
    virtual double get_x_jump() const = 0;

    /// @brief Начальная плотность
    virtual double density(const Vector3d &r) const = 0;

    /// @brief Начальная скорость
    virtual Vector3d velocity(const Vector3d &r) const = 0;

    /// @brief Начальное давление
    virtual double pressure(const Vector3d &r) const = 0;

    /// @brief Начальная внутренняя энергия
    virtual double energy(const Vector3d &r) const = 0;

    /// @brief Уравнение состояния
    /// (для одноматериальных тестов не зависит от r)
    virtual Eos::Ptr get_eos(const Vector3d &r) const = 0;

    virtual ~Test1D() = default;
};

} // namespace zephyr::phys