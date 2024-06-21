#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/eos.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;

/// @class Абстрактный класс одномерного теста
class Test1D {
public:
    /// @brief Получить название теста
    virtual std::string get_name() const = 0;

    /// @brief Левая граница области
    virtual double xmin() const = 0;

    /// @brief Правая граница области
    virtual double xmax() const = 0;

    /// @brief Конечный момент времени
    virtual double max_time() const = 0;

    /// @brief Получить используемый УрС
    virtual const Eos& get_eos() const = 0;

    /// @brief Получить положение разрыва
    virtual double get_x_jump() const = 0;


    /// @brief Начальная плотность
    virtual double density(double x) const = 0;

    /// @brief Начальная скорость
    virtual Vector3d velocity(double x) const = 0;

    /// @brief Начальное давление
    virtual double pressure(double x) const = 0;

    /// @brief Начальная внутренняя энергия
    virtual double energy(double x) const = 0;

    /// @brief Версия для многоматериальных
    virtual const Eos& get_eos(double x) const {
        return get_eos();
    }


    /// @brief Начальная плотность
    virtual double density(const Vector3d &r) const {
        return density(r.x());
    }

    /// @brief Начальная скорость
    virtual Vector3d velocity(const Vector3d &r) const {
        return velocity(r.x());
    }

    /// @brief Начальное давление
    virtual double pressure(const Vector3d &r) const {
        return pressure(r.x());
    }

    /// @brief Начальная внутренняя энергия
    virtual double energy(const Vector3d &r) const {
        return energy(r.x());
    }

    /// @brief Версия для многоматериальных
    virtual const Eos& get_eos(const Vector3d &r) const {
        return get_eos(r.x());
    }

    virtual ~Test1D() = default;
};

}