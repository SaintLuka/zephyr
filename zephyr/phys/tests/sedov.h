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
    double V_init; ///<

    /// @brief Конструктор
    SedovBlast(double gamma=1.4) : eos(gamma) {
        
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
        return r.norm() < x_jump ? r_init : r_amb;
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &r) const { 
        return r.norm() < x_jump ? v_init : v_amb;
    }

    /// @brief Начальное давление
    double pressure(const Vector3d &r) const { 
        return r.norm() < x_jump ? p_init : p_amb;
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const { 
        return eos.energy_rp(density(r), pressure(r));
    }

    ~SedovBlast() = default;
};

}
}