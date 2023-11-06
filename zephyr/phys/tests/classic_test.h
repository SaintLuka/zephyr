#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;

/// @class Абстрактный класс для классического одномерного теста
class ClassicTest {
public:
    /// @brief Получить название теста
    [[nodiscard]] virtual std::string get_name() const = 0;

    /// @brief Левая граница области
    [[nodiscard]] virtual double xmin() const = 0;

    /// @brief Правая граница области
    [[nodiscard]] virtual double xmax() const = 0;

    /// @brief Конечный момент времени
    [[nodiscard]] virtual double max_time() const = 0;

    ///@brief Получить используемый УрС
    [[nodiscard]] virtual const Eos& get_eos() const = 0;

    ///@brief Получить положение разрыва
    [[nodiscard]] virtual double get_x_jump() const = 0;

    /// @brief Начальная плотность
    [[nodiscard]] virtual double density(const double &x) const = 0;

    /// @brief Начальная скорость
    [[nodiscard]] virtual Vector3d velocity(const double &x) const = 0;

    /// @brief Начальное давление
    [[nodiscard]] virtual double pressure(const double &x) const = 0;

    /// @brief Начальная внутренняя энергия
    [[nodiscard]] virtual double energy(const double &x) const = 0;


    /// @brief Начальная плотность
    [[nodiscard]] virtual double density(const Vector3d &r) const = 0;

    /// @brief Начальная скорость
    [[nodiscard]] virtual Vector3d velocity(const Vector3d &r) const = 0;

    /// @brief Начальное давление
    [[nodiscard]] virtual double pressure(const Vector3d &r) const = 0;

    /// @brief Начальная внутренняя энергия
    [[nodiscard]] virtual double energy(const Vector3d &r) const = 0;

    virtual ~ClassicTest() = default;
};

}