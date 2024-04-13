#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr {
namespace phys {

using zephyr::geom::Vector3d;

class ShockInABox {
public:

    IdealGas eos;   ///< Используемый УрС
    double r_jump; //m ///< Положение разрыва
    double finish = 1.25 * 1e-3;  ///< Конечный момент времени
    double p0, p1; //atm ///<Давление
    double r0, r1; // kg/m^3 ///< Плотность
    Vector3d v0, v1;  ///< Скорость

    /// @brief Конструктор
    // Ms - число Маха
    ShockInABox(double Ms = 4.5) : eos(1.4) {    

        // 1.4979066.pdf
        // Initial conditions for the shock-in-a-box case
        // On the propagation and multiple reflections of a blast wave travelling through a dusty gas in a closed box

        /*
        With the pre-shock state corresponding to the properties
        of the undisturbed gas (po, ρo), we determine the post-shock
        thermodynamic state using the standard Rankine–Hugoniot
        conditions for a travelling shock 
        */
        p0 = 1.0;
        r0 = 1.225;
        v0 = { 0, 0, 0 };

        r1 = r0 * ( eos.gamma + 1 ) * Ms * Ms / ( 2 + ( eos.gamma - 1 ) * Ms * Ms );
        p1 = p0 * ( 1 + 2 * eos.gamma / ( eos.gamma + 1 ) * ( Ms * Ms  - 1) );
        double vr = Ms * std::sqrt( eos.gamma * p0 * r0 ) * ( 1 - r0 / p0 ); // V1 is the radial velocity imparted by the shock wave to the fluid at r = rs.
        v1 = { vr / std::sqrt( 2 ), vr / std::sqrt( 2 ), 0 };

        // We initially place the shock at a radial distance rs = 0.1L from the origin of the blast
        r_jump = 0.5;
        finish = 0.4;
    }

    std::string get_name() const { return "ShockInABox";}

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
    double get_r_jump() const { return r_jump; }

    /// @brief Начальная плотность
    double density(const Vector3d &rs) const 
    { 
        Vector3d rs2 = rs - Vector3d{0.5, 0.5, 0};
        return rs2.norm() > r_jump ? r0 : r1; 
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &rs) const 
    { 
        Vector3d rs2 = rs - Vector3d{0.5, 0.5, 0};
        return rs2.norm() > r_jump ? v0 : v1;
    }

    /// @brief Начальное давление
    double pressure(const Vector3d &rs) const 
    { 
        Vector3d rs2 = rs - Vector3d{0.5, 0.5, 0};
        return rs2.norm() > r_jump ? p0 : p1;
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &rs) const { 
        return eos.energy_rp(density(rs), pressure(rs));
    }

    ~ShockInABox() = default;
};

}
}