#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr {
namespace phys {

using zephyr::geom::Vector3d;

/// @class Абстрактный класс для классического одномерного теста
class ClassicTest {
public:
    /// @brief Получить название теста
    virtual std::string get_name() const = 0;

    /// @brief Левая граница области
    virtual double xmin() const = 0;

    /// @brief Правая граница области
    virtual double xmax() const = 0;

    /// @brief Конечный момент времени
    virtual double max_time() const = 0;

    ///@brief Получить используемый УрС
    virtual const Eos* get_eos() const = 0;

    ///@brief Получить положение разрыва
    virtual double get_x_jump() const = 0;

    /// @brief Начальная плотность
    virtual double density(const double &x) const = 0;

    /// @brief Начальная скорость
    virtual Vector3d velocity(const double &x) const = 0;

    /// @brief Начальное давление
    virtual double pressure(const double &x) const = 0;

    /// @brief Начальная внутренняя энергия
    virtual double energy(const double &x) const = 0;


    /// @brief Начальная плотность
    virtual double density(const Vector3d &r) const = 0;

    /// @brief Начальная скорость
    virtual Vector3d velocity(const Vector3d &r) const = 0;

    /// @brief Начальное давление
    virtual double pressure(const Vector3d &r) const = 0;

    /// @brief Начальная внутренняя энергия
    virtual double energy(const Vector3d &r) const = 0;

    virtual ~ClassicTest() = default;
};

}
}