#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr {
namespace phys {

using zephyr::geom::Vector3d;

class RiemannTest2D {
public:

    IdealGas m_eos;   ///< Используемый УрС
    double x_jump;  ///< Положение разрыва
    double y_jump;

    double finish;  ///< Конечный момент времени
    
    double r1, r2, r3, r4;  ///< Плотность
    double u1, u2, u3, u4;  ///< Скорость x
    double v1, v2, v3, v4;  ///< Скорость y
    double p1, p2, p3, p4;  ///< Давление
    double e1, e2, e3, e4;  ///< Внутренняя энергия


    /// @brief Конструктор
    explicit RiemannTest2D (int testCase) : m_eos(1.4) {
        switch (testCase)
        {
        case 1: {
            r1 = 1; r2 = 0.5197; r3 = 0.1072; r4 = 0.2579;
            u1 = 0; u2= -0.7259; u3 = -0.7259; u4 = 0;
            v1 = 0; v2 = 0; v3 = -1.4045; v4 = -1.4045;
            p1 = 1; p2 = 0.4; p3 = 0.0439; p4 = 0.15;
            break;
        }
        case 2:
            break;
        case 3:
            break;

        
        default:
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

    /// @brief Конечный момент времени
    double max_time() const { return finish; }

    ///@brief Получить используемый УрС
    const Eos& get_eos() const { return m_eos; }

    ///@brief Получить положение разрыва
    double get_x_jump() const { return x_jump; }

    ///@brief Получить положение разрыва
    double get_y_jump() const { return y_jump; }

    /// @brief Начальная плотность
    double density(const Vector3d &vec) const { 
        if (vec.x() > 0.5 && vec.y() > 0.5)
            return r1;
        if (vec.x() < 0.5 && vec.y() > 0.5)
            return r2;
        if (vec.x() < 0.5 && vec.y() < 0.5)
            return r3;
        if (vec.x() > 0.5 && vec.y() < 0.5)
            return r4;
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &vec) const { 
        if (vec.x() > 0.5 && vec.y() > 0.5)
            return {u1, v1, 0.0};
        if (vec.x() < 0.5 && vec.y() > 0.5)
            return {u2, v2, 0.0};
        if (vec.x() < 0.5 && vec.y() < 0.5)
            return {u3, v3, 0.0};
        if (vec.x() > 0.5 && vec.y() < 0.5)
            return {u4, v4, 0.0};
    }

    /// @brief Начальное давление
    double pressure(const Vector3d &vec) const { 
        if (vec.x() > 0.5 && vec.y() > 0.5)
            return p1;
        if (vec.x() < 0.5 && vec.y() > 0.5)
            return p2;
        if (vec.x() < 0.5 && vec.y() < 0.5)
            return p3;
        if (vec.x() > 0.5 && vec.y() < 0.5)
            return p4;
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const { 
        return m_eos.energy_rp(density(r), pressure(r));
    }

    ~RiemannTest2D() = default;
};

}
}