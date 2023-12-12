#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr {
namespace phys {

using zephyr::geom::Vector3d;

class BlastWave {
public:

    IdealGas eos;   ///< Используемый УрС
    double r_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rS, ra;  ///< Плотность
    Vector3d vS, va;  ///< Скорость
    double pS, pa;  ///< Давление
    double eS, ea;  ///< Внутренняя энергия
    int type = 1;
    Vector3d center = {0.5, 0.5, 0};


    /// @brief Конструктор
    BlastWave(int type) : eos(1.4) {
        
        rS = 0.125;
        vS = Vector3d{0.0, 0.0, 0.0};
        pS = 0.1;

        ra = 1.0;
        va = Vector3d{0.0, 0.0, 0.0};
        pa = 1.0;
        
        r_jump = 0.2;
        finish = 1.0;
        type = type;;
    }

    std::string get_name() const { return "BlastWave";}

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
    double density(const Vector3d &x) const 
    { 
        if (type == 1) {
            double r = x.norm() / r_jump;
            return r < 1.0 ? rS : ra; 
        }
        if (type == 2) {
            double Mach = 0;
        }
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &x) const 
    { 
        if (type == 1) {
            double r = x.norm() / r_jump;
            return r < 1.0 ? vS : va;
        }
        if (type == 2) {
        }
    }

    /// @brief Начальное давление
    double pressure(const Vector3d &x) const 
    { 
        if (type == 1) {
            double r = x.norm() / r_jump;
            return r < 1.0 ? pS : pa; 
        }
        if (type == 2) {

        }
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &x) const { 
        return eos.energy_rp(density(x), pressure(x));
    }

    ~BlastWave() = default;
};

}
}