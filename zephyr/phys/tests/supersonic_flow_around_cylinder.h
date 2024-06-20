#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/tests/classic_test.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;

/// @class Классический одномерный тест Сода
class SuperSonicFlowAroundCylinder : public ClassicTest {
public:
    IdealGas eos;   ///< Используемый УрС
    double finish;  ///< Конечный момент времени
    double rho;  ///< Плотность
    double u, v;  ///< Скорость
    double p;  ///< Давление
    double e;  ///< Внутренняя энергия
    double x_jump = 0;

    // Bow shock resulting from a supersonic flow around a stationary cylinder 
    // tends to suffer from the carbuncle phenomenon
    // Peery and Imlay, 1988

    /// @brief Конструктор
    SuperSonicFlowAroundCylinder(double Ms=3) : eos(1.4) {
        // @formatter:off
        rho = 1;
        u = sqrt(1.4) * Ms;
        v = 0;
        p = 1;

        finish = 10;
        // @formatter:on
    }

    std::string get_name() const { return "Supersonic flow around cylinder";}

    /// @brief Левая граница области
    double xmin() const { return 0.0; }

    /// @brief Правая граница области
    double xmax() const { return 1.0; }

    /// @brief 
    double ymin() const { return 0.0; }

    /// @brief  
    double ymax() const {return 0.8; }

    /// @brief Конечный момент времени
    double max_time() const { return finish; }

    /// @brief  
    double get_x_jump() const { return x_jump; }

    ///@brief Получить используемый УрС
    const Eos& get_eos() const { return eos; }

    /// @brief Начальная плотность
    double density(const double &x) const { return rho; }
    double density(const Vector3d &r) const { return rho; }

    /// @brief Начальная внутренняя энергия
    double energy(const double &x) const { return eos.energy_rp(rho, p); }
    double energy(const Vector3d &r) const { return eos.energy_rp(rho, p); }
    
    /// @brief Начальное давление
    double pressure(const double &x) const { return p; }
    double pressure(const Vector3d &r) const { return p; }

    /// @brief Начальная скорость
    Vector3d velocity(const double &x) const { return Vector3d(u, v, 0); }
    Vector3d velocity(const Vector3d &r) const { return Vector3d(u, v, 0); }


    ~SuperSonicFlowAroundCylinder() override = default;

};

}

