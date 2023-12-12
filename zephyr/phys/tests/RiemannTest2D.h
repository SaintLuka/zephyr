#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr {
namespace phys {

using zephyr::geom::Vector3d;

class RiemannTest2D {
public:

    IdealGas eos;   ///< Используемый УрС
    double x_jump;  ///< Положение разрыва
    double y_jump;

    double finish;  ///< Конечный момент времени
    
    double r1, r2, r3, r4;  ///< Плотность
    double u1, u2, u3, u4;  ///< Скорость x
    double v1, v2, v3, v4;  ///< Скорость y
    double p1, p2, p3, p4;  ///< Давление
    double e1, e2, e3, e4;  ///< Внутренняя энергия


    /// @brief Конструктор
    explicit RiemannTest2D (int testCase) : eos(1.4) {
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
            u1 = -0.75; u2 = 0.75;  u3 = -0.75; u4 = -0.75;
            v1 = -0.5;  v2 = 0.5;   v3 = 0.5;   v4 = -0.5;
            p1 = 1;     p2 = 1;     p3 = 1;     p4 = 1;
            break; 
        }

        y_jump = 0.5;
        x_jump = 0.5;
        finish = 1.0;
    };

    std::string get_name() const { return "2D Riemann Test";}

    /// @brief Левая граница области
    double xmin() const { return 0.0; }

    /// @brief Правая граница области
    double xmax() const { return 1.0; }

    /// @brief Левая граница области
    double ymin() const { return 0.0; }

    /// @brief Правая граница области
    double ymax() const { return 1.0; }

    /// @brief Конечный момент времени
    double max_time() const { return finish; }

    ///@brief Получить используемый УрС
    const Eos& get_eos() const { return eos; }

    ///@brief Получить положение разрыва
    double get_x_jump() const { return x_jump; }

    ///@brief Получить положение разрыва
    double get_y_jump() const { return y_jump; }

    /// @brief Начальная плотность
    double density(const Vector3d &vec) const { 
        if (vec.x() > x_jump && vec.y() > y_jump)
            return r1;
        if (vec.x() < x_jump && vec.y() > y_jump)
            return r2;
        if (vec.x() < x_jump && vec.y() < y_jump)
            return r3;
        if (vec.x() > x_jump && vec.y() < y_jump)
            return r4;
    };

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &vec) const { 
        if (vec.x() > x_jump && vec.y() > y_jump)
            return {u1, v1, 0.0};
        if (vec.x() < x_jump && vec.y() > y_jump)
            return {u2, v2, 0.0};
        if (vec.x() < x_jump && vec.y() < y_jump)
            return {u3, v3, 0.0};
        if (vec.x() > x_jump && vec.y() < y_jump)
            return {u4, v4, 0.0};
    };

    /// @brief Начальное давление
    double pressure(const Vector3d &vec) const { 
        if (vec.x() > x_jump && vec.y() > y_jump)
            return p1;
        if (vec.x() < x_jump && vec.y() > y_jump)
            return p2;
        if (vec.x() < x_jump && vec.y() < y_jump)
            return p3;
        if (vec.x() > x_jump && vec.y() < y_jump)
            return p4;
    };

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const { 
        return eos.energy_rp(density(r), pressure(r));
    };

    ~RiemannTest2D() = default;
};

}
}