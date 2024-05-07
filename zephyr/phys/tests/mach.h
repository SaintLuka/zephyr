#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/tests/classic_test.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;

/// @class Классический одномерный тест Сода
class Mach : public ClassicTest {
public:
    IdealGas eos;   ///< Используемый УрС
    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия

    /// Double Mach  Reflection of a Strong Shock
    /// Woodward and Colella

    /// @brief Конструктор
    Mach(double Ms=1.2) : eos(1.4) {
        // @formatter:off
        pR = 1.01325;
        rR = 1.225;
        uR = 0;
        
        pL = pR * (2 * eos.gamma * Ms * Ms  - eos.gamma + 1) / (eos.gamma + 1) ;
        rL = rR * ( eos.gamma + 1 ) * Ms * Ms / ( 2 + ( eos.gamma - 1 ) * Ms * Ms );
        uL = 2 / Ms * std::sqrt( eos.gamma * pR / rR ) * ( Ms * Ms - 1 ) / (eos.gamma + 1 ); 

        eL = eos.energy_rp(rL, pL);
        eR = eos.energy_rp(rR, pR);

        x_jump = 0.1;
        finish = 0.5;
        // @formatter:on
    }

    /// @brief Симметрично отразить начальные условия
    void inverse() {
        std::swap(rL, rR);
        std::swap(uL, uR);
        std::swap(pL, pR);
        std::swap(eL, eR);
        uL *= -1.0;
        uR *= -1.0;
    }

    std::string get_name() const { return "Mach";}

    /// @brief Левая граница области
    double xmin() const { return 0.0; }

    /// @brief Правая граница области
    double xmax() const { return 0.5; }

    /// @brief Конечный момент времени
    double max_time() const { return finish; }

    ///@brief Получить используемый УрС
    const Eos& get_eos() const { return eos; }

    ///@brief Получить положение разрыва
    double get_x_jump() const { return x_jump; }

    /// @brief Начальная плотность
    double density(const double &x) const { return x < x_jump ? rL : rR; }

    /// @brief Начальная скорость
    Vector3d velocity(const double &x) const { return { x < x_jump ? uL : uR, 0.0, 0.0}; }

    /// @brief Начальное давление
    double pressure(const double &x) const { return x < x_jump ? pL : pR; }

    /// @brief Начальная внутренняя энергия
    double energy(const double &x) const { return x < x_jump ? eL : eR; }

    /// @brief Начальная плотность
    double density(const Vector3d &r) const { return density(r.x()); }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &r) const { return velocity(r.x()); }

    /// @brief Начальное давление
    double pressure(const Vector3d &r) const { return pressure(r.x()); }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const { return energy(r.x()); }

    ~Mach() override = default;

};

}

