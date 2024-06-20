#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr {
namespace phys {

using zephyr::geom::Vector3d;

class SedovBlast {
public:

    IdealGas eos;   ///< Используемый УрС
    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double r_init, r_amb;  ///< Плотность
    Vector3d v_init, v_amb;  ///< Скорость
    double p_init, p_amb;  ///< Давление
    double E; ///< Начальная энергия
    double V_init; ///< Начальный объем
    double Ms; ///< Число Маха
    int constructor_type;


    /// @brief Конструктор
    SedovBlast(double gamma=1.4) : eos(gamma), constructor_type(1) {
        // The Sedov (or Sedov-Taylor) blast wave is a standard hydrodynamics test problem
        // A large amount of energy is placed into a very small volume, driving a spherical 
        // (or cylindrical in 2-d Cartesian coordinates) blast wave
        
        E = 1; // у.е.
        x_jump = 0.02;
        V_init = M_PI * x_jump * x_jump;
        
        p_amb = 0.06;
        p_init = (eos.gamma - 1) * E / V_init;
        
        r_amb = 1.0;
        r_init = 1.0;

        v_amb = {0,0,0};
        v_init = {0,0,0};
        
        finish = 1.0;
    }

    /// @brief  задать число Маха и начальный радиус ударной волны
    SedovBlast(double Ms, double x_jump, double gamma=1.4) : eos(gamma), x_jump(x_jump), constructor_type(2) {
        r_amb = 1.225;
        p_amb = 1;
        v_amb = {0,0,0};

        r_init = r_amb * (eos.gamma + 1) * Ms * Ms / (2 + (eos.gamma - 1) * Ms * Ms);
        p_init = p_amb * (1 + 2 * eos.gamma * (Ms * Ms - 1) / (eos.gamma + 1));
        double radV = Ms * sqrt(eos.gamma * p_amb / r_amb) * (1 - r_amb / r_init);
        v_init = {radV / sqrt(2), radV / sqrt(2), 0};

        finish = 1.0;
    }

    std::string get_name() const { return "SedovBlast";}

    /// @brief Левая граница области
    double xmin() const { return -0.5; }

    /// @brief Правая граница области
    double xmax() const { return 0.5; }

    /// @brief Левая граница области
    double ymin() const { return -0.5; }

    /// @brief Правая граница области
    double ymax() const { return 0.5; }

    /// @brief Конечный момент времени
    double max_time() const { return finish; }

    ///@brief Получить используемый УрС
    const Eos& get_eos() const { return eos; }

    ///@brief Получить положение разрыва
    double get_x_jump() const { return x_jump; }

    /// @brief Начальная плотность
    double density(const Vector3d &r) const { 
        switch (constructor_type)
        {
        case 1:
            return r.norm() < x_jump ? r_init : r_amb;
        case 2:
            return r.norm() < x_jump ? std::pow((r.norm() / x_jump), 3 / (eos.gamma - 1)) : r_amb;
        default:
            break;
        }    
    }

    /// @brief Начальное давление
    double pressure(const Vector3d &r) const { 
        switch (constructor_type)
        {
        case 1:
            return r.norm() < x_jump ? p_init : p_amb;
        case 2:
            return r.norm() < x_jump ? p_init * std::pow(3, std::pow(r.norm() / x_jump, 5)) : p_amb;
        default:
            break;
        }
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &r) const { 
        switch (constructor_type)
        {
        case 1:
            return r.norm() < x_jump ? v_init : v_amb;
        case 2:
            return r.norm() < x_jump ? v_init * (r.norm() / x_jump) : v_amb;
        default:
            break;
        }
    }


    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const { 
        return eos.energy_rp(density(r), pressure(r));
    }

    ~SedovBlast() = default;
};

}
}