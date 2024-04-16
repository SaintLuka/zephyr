#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr {
namespace phys {

using zephyr::geom::Vector3d;

class BlastWave {
public:

    IdealGas eos;   ///< Используемый УрС
    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rS, ra;  ///< Плотность
    Vector3d vS, va;  ///< Скорость
    double pS, pa;  ///< Давление
    double eS, ea;  ///< Внутренняя энергия
    double V1; ///< Радиальная скорость


    /// @brief Конструктор
    BlastWave(double Ms) : eos(1.4) {
        ra = 1.225;
        va = Vector3d{0.0, 0.0, 0.0};
        pa = 1.0;

        pS = pa * ( 1 + 2 * eos.gamma / ( eos.gamma + 1 ) * ( Ms * Ms  - 1) );
        rS = ra * ( eos.gamma + 1 ) * Ms * Ms / ( 2 + ( eos.gamma - 1 ) * Ms * Ms );
        V1 = Ms * std::sqrt( eos.gamma * pa / ra ) * ( 1 - ra / rS );

        x_jump = 0.1;
        finish = 1.0;
    }

    std::string get_name() const { return "BlastWave";}

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
    double density(const Vector3d &x) const 
    { 
        double ratio = x.norm() / x_jump;
        if (ratio < 1)
            return rS * 0.5 * (std::exp(std::pow(ratio, 5) - 1));
        return ra; 
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &x) const 
    { 
        double ratio = x.norm() / x_jump;
        if (ratio < 1) 
            return Vector3d{V1 * x.x() / x.norm(), V1 * x.y() / x.norm(), 0};
        return va; 
    }

    /// @brief Начальное давление
    double pressure(const Vector3d &x) const 
    { 
        double ratio = x.norm() / x_jump;
        if (ratio < 1)
            return pS * 0.3 * (std::exp(std::pow(ratio, 5)));
        return pa;
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &x) const { 
        return eos.energy_rp(density(x), pressure(x));
    }

    ~BlastWave() = default;
};

}
}