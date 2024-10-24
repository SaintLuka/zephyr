#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/geom/primitives/boundary.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;
using zephyr::geom::Boundary;
using zephyr::geom::generator::Rectangle;

/// @brief Двумерный распад разрыва.
/// Ссылка на статью???
class Test2D {
public:

    IdealGas::Ptr eos;  ///< Используемый УрС
    double x_jump;      ///< Положение разрыва
    double y_jump;      ///< Положение разрыва

    double finish;      ///< Конечный момент времени
    
    double r1, r2, r3, r4;  ///< Плотность
    double u1, u2, u3, u4;  ///< Скорость x
    double v1, v2, v3, v4;  ///< Скорость y
    double p1, p2, p3, p4;  ///< Давление
    double e1, e2, e3, e4;  ///< Внутренняя энергия

    Rectangle generator;    ///< Сеточный генератор


    /// @brief Конструктор
    explicit Test2D (int testCase) {
        eos = IdealGas::create(1.4);

        switch (testCase)
        {
        case 1: {
            r1 = 1; r2 = 0.5197;  r3 = 0.1072;  r4 = 0.2579;
            u1 = 0; u2 = -0.7259; u3 = -0.7259; u4 = 0;
            v1 = 0; v2 = 0;       v3 = -1.4045; v4 = -1.4045;
            p1 = 1; p2 = 0.4;     p3 = 0.0439;  p4 = 0.15;
            break;
        }
        case 2:
            r1 = 1; r2 = 0.4;     r3 = 1;       r4 = 0.5197;
            u1 = 0; u2 = -0.7259; u3 = -0.7259; u4 = 0;
            v1 = 0; v2 = 0;       v3 = -0.7259; v4 = -0.7259;
            p1 = 1; p2 = 0.5197;  p3 = 1;       p4 = 0.4;
            break;
        case 3:
            r1 = 1.5; r2 = 0.5323; r3 = 0.138; r4 = 0.5323;
            u1 = 0;   u2= 1.206;   u3 = 1.206; u4 = 0;
            v1 = 0;   v2 = 0;      v3 = 1.206; v4 = 1.206;
            p1 = 1.5; p2 = 0.3;    p3 = 0.029; p4 = 0.3;
            break;
        case 4:
            r1 = 1.1; r2 = 0.5065; r3 = 1.1;    r4 = 0.5065;
            u1 = 0;   u2= 0.8939;  u3 = 0.8939; u4 = 0;
            v1 = 0;   v2 = 0;      v3 = 0.8939; v4 = 0.8939;
            p1 = 1.1; p2 = 0.35;   p3 = 1.1;    p4 = 0.35;
            break;
        case 5:
            r1 = 1;     r2 = 2;     r3 = 1;    r4 = 3;
            u1 = -0.75; u2 = -0.75; u3 = 0.75; u4 = 0.75;
            v1 = -0.5;  v2 = 0;     v3 = 0.5;  v4 = -0.5;
            p1 = 1;     p2 = 1;     p3 = 1;    p4 = 1;
            break; 
        case 6:
            r1 = 1;     r2 = 2;     r3 = 1;     r4 = 3;
            u1 = 0.75;  u2 = 0.75;  u3 = -0.75; u4 = -0.75;
            v1 = -0.5;  v2 = 0.5;   v3 = 0.5;   v4 = -0.5;
            p1 = 1;     p2 = 1;     p3 = 1;     p4 = 1;
            finish = 0.3;
            break; 
        case 7:
            r1 = 0.5197; r2 = 1.0; r3 = 0.8; r4 = 1.0;
            p1 = 0.4; p2 = 1.0; p3 = 1.0; p4 = 1.0;
            u1 = 0.1; u2 = -0.6259; u3 = 0.1; u4 = 0.1;
            v1 = 0.1; v2 = 0.1; v3 = 0.1; v4 = -0.6259;
            finish = 0.25;
            break;
        case 8:
            r1 = 1.0; r2 = 2.0; r3 = 1.039; r4 = 0.5197;
            p1 = 1.0; p2 = 1.0; p3 = 0.4; p4 = 0.4;
            u1 = 0.0; u2 = 0.0; u3 = 0.0; u4 = 0.0;
            v1 = 0.3; v2 = -0.3; v3 = -0.8133; v4 = -0.259;
            break;
        default:
            throw std::runtime_error("Unknown Test");
        }

        y_jump = 0.5;
        x_jump = 0.5;

        generator = Rectangle(0.0, 1.0, 0.0, 1.0);
        generator.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
                                  .bottom=Boundary::ZOE, .top=Boundary::ZOE});
    };

    /// @brief Название теста
    std::string name() const { return "2D Riemann Test";}

    /// @brief Конечный момент времени
    double max_time() const { return finish; }

    /// @brief Начальная плотность
    double density(const Vector3d &vec) const {
        if (vec.x() >= x_jump && vec.y() >= y_jump)
            return r1;
        if (vec.x() <= x_jump && vec.y() >= y_jump)
            return r2;
        if (vec.x() <= x_jump && vec.y() <= y_jump)
            return r3;
        if (vec.x() >= x_jump && vec.y() <= y_jump)
            return r4;
    };

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &vec) const { 
        if (vec.x() >= x_jump && vec.y() >= y_jump)
            return {u1, v1, 0.0};
        if (vec.x() <= x_jump && vec.y() >= y_jump)
            return {u2, v2, 0.0};
        if (vec.x() <= x_jump && vec.y() <= y_jump)
            return {u3, v3, 0.0};
        if (vec.x() >= x_jump && vec.y() <= y_jump)
            return {u4, v4, 0.0};
    };

    /// @brief Начальное давление
    double pressure(const Vector3d &vec) const { 
        if (vec.x() >= x_jump && vec.y() >= y_jump)
            return p1;
        if (vec.x() <= x_jump && vec.y() >= y_jump)
            return p2;
        if (vec.x() <= x_jump && vec.y() <= y_jump)
            return p3;
        if (vec.x() >= x_jump && vec.y() <= y_jump)
            return p4;
    };

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const { 
        return eos->energy_rP(density(r), pressure(r));
    };

    Eos::Ptr get_eos(const Vector3d& r) const {
        return eos;
    }

    ~Test2D() = default;
};

} // namespace zephyr::phys