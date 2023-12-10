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
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия


    /// @brief Конструктор
    BlastWave() : eos(1.4) {
        rL = 0.125;
        uL = 0.0;
        pL = 0.1;
        rR = 1.0;
        uR = 0.0;
        pR = 1.0;
        x_jump = 0.4;
        finish = 1.0;
    }

    std::string get_name() const { return "BlastWave";}

    /// @brief Левая граница области
    double xmin() const { return 0.0; }

    /// @brief Правая граница области
    double xmax() const { return 2.0; }

    /// @brief Конечный момент времени
    double max_time() const { return finish; }

    ///@brief Получить используемый УрС
    const Eos& get_eos() const { return eos; }

    ///@brief Получить положение разрыва
    double get_x_jump() const { return x_jump; }

    /// @brief Начальная плотность
    double density(const Vector3d &x) const 
    { 
        double r = x.norm() / x_jump;
        return r < 1.0 ? rL : rR; 
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &x) const { 
        return {0.0, 0.0, 0.0}; 
    }

    /// @brief Начальное давление
    double pressure(const Vector3d &x) const 
    { 
        double r = x.norm() / x_jump;
        return r < 1.0 ? pL : pR; 
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &x) const { 
        return eos.energy_rp(density(x), pressure(x));
    }

    ~BlastWave() = default;
};

}
}