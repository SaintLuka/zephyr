#if 0
#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/tests/test_1D.h>

#include <zephyr/geom/generator/rectangle.h>

namespace zephyr::phys {

using namespace zephyr::geom;
using zephyr::geom::generator::Rectangle;


/// @brief Косая ударная волна
class SkewShockWave {
public:
    ShockWave sw;
    double angle;         ///< Угол поворота
    Rectangle generator;  ///< Сеточный генератор

    /// @brief Конструктор
    SkewShockWave(double Ms=6, double angle = M_PI/6, double finish=0.2)
            : sw(Ms, 0.0, finish), angle(angle) {

        generator = Rectangle(-0.2, 2.3, 0.0, 1.0);
        generator.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
                                  .bottom=Boundary::WALL, .top=Boundary::ZOE});

    }

    std::string name() const { return "Skew Shock Wave";}

    /// @brief Левая граница области
    double xmin() const { return -0.2; }

    /// @brief Правая граница области
    double xmax() const { return 2.3; }

    /// @brief Конечный момент времени
    double max_time() const { return sw.finish; }

    /// @brief Вращение координат
    Vector3d rotate(const Vector3d& r, double phi) const {
        Vector3d res = {
                +std::cos(phi) * r.x() - std::sin(phi) * r.y(),
                +std::sin(phi) * r.x() + std::cos(phi) * r.y(),
                0.0
        };
        return res;
    }


    /// @brief Начальная плотность
    double density(const Vector3d& r) const {
        return sw.density(rotate(r, angle));
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d& r) const {
        return rotate(sw.velocity(rotate(r, angle)), -angle);
    }

    /// @brief Начальное давление
    double pressure(const Vector3d& r) const {
        return sw.pressure(rotate(r, angle));
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d& r) const {
        return sw.energy(rotate(r, angle));
    }

    /// @brief Уравнение состояния
    Eos::Ptr get_eos(const Vector3d& r) const {
        return sw.get_eos(r);
    }
};

} // namespace zephyr::phys

#endif