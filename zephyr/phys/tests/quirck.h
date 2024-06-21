#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/tests/test_1D.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;

class QuirckTest : public Test1D {
public:
    IdealGas eos;   ///< Используемый УрС
    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия

    /// @brief Конструктор
    QuirckTest(double Ms=10, double x_jump=0.1, double finish=0.1) : eos(1.4), x_jump(x_jump), finish(finish) {
        // @formatter:off
        pR = 1.0;
        rR = 1.0;
        uR = 0;

        /// Double Mach  Reflection of a Strong Shock
        /// Woodward and Colella ????
        
        pL = pR * (2 * eos.gamma * Ms * Ms  - eos.gamma + 1) / (eos.gamma + 1) ;
        rL = rR * ( eos.gamma + 1 ) * Ms * Ms / ( 2 + ( eos.gamma - 1 ) * Ms * Ms );
        uL = 2 / Ms * std::sqrt( eos.gamma * pR / rR ) * ( Ms * Ms - 1 ) / (eos.gamma + 1 ); 

        eL = eos.energy_rp(rL, pL);
        eR = eos.energy_rp(rR, pR);
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

    std::string get_name() const final { return "QuirckTest";}

    /// @brief Левая граница области
    double xmin() const final { return 0.0; }

    /// @brief Правая граница области
    double xmax() const final { return 1.5; }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    ///@brief Получить используемый УрС
    const Eos& get_eos() const final { return eos; }

    ///@brief Получить положение разрыва
    double get_x_jump() const final { return x_jump; }


    /// @brief Начальная плотность
    double density(double x) const final {
        return x < x_jump ? rL : rR;
    }

    /// @brief Начальная скорость
    Vector3d velocity(double x) const final {
        return { x < x_jump ? uL : uR, 0, 0};
    }

    /// @brief Начальное давление
    double pressure(double x) const final {
        return x < x_jump ? pL : pR;
    }

    /// @brief Начальная внутренняя энергия
    double energy(double x) const final {
        return x < x_jump ? eL : eR;
    }


    /// @brief Начальная плотность
    double density(const Vector3d& r) const final {
        return density(r.x());
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d& r) const final {
        return velocity(r.x());
    }

    /// @brief Начальное давление
    double pressure(const Vector3d& r) const final {
        return pressure(r.x());
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d& r) const final {
        return energy(r.x());
    }

    ~QuirckTest() final = default;

};

}

