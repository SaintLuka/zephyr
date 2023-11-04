#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/tests/classic_test.h>

namespace zephyr { 
namespace phys {

using zephyr::geom::Vector3d;

/// @brief Тест Шу-Ошера 
/// C.-W. Shu and S. Osher. Efficient Implementation of Essentially 
/// Non-oscillatory Shock-Capturing Schemes, II (1988)
class ShuOsherTest : public ClassicTest {
public:
    IdealGas eos;   ///< Используемый УрС
    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия
    double epsilon; ///< Эпсилон

    /// @brief Начальные данные
    ShuOsherTest() : eos(1.4) {
        rL  = 27.0 / 7.0; rR  = 1.0;
        pL = 31.0 / 3.0; pR = 1.0;
        uL = 4.0 * std::sqrt(35.0) / 9.0; uR = 0.0;
        eL = eos.energy_rp(rL, pL);
        eR = eos.energy_rp(rR, pR);

        x_jump = -4.0;
        finish = 2.0;
        epsilon = 0.2;
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

    std::string get_name() const { return "ShuOsherTest";}

    /// @brief Левая граница области
    double xmin() const { return -5.0; }

    /// @brief Правая граница области
    double xmax() const { return +5.0; }

    /// @brief Конечный момент времени
    double max_time() const { return finish; }

    /// @brief Получить используемый УрС
    const Eos& get_eos() const { return eos; }

    ///@brief Получить положение разрыва
    double get_x_jump() const { return x_jump; }

    /// @brief Начальная плотность
    double density(const double &x) const { 
        return x < x_jump ? rL : rR + epsilon * std::sin(5.0 * x); 
    }

    /// @brief Начальная скорость
    Vector3d velocity(const double &x) const { 
        return {x < x_jump ? uL : uR, 0.0, 0.0}; 
    }
    
    /// @brief Начальное давление
    double pressure(const double &x) const { 
        return x < x_jump ? pL : pR; 
    }

    /// @brief Начальная внутренняя энергия
    double energy(const double &x) const { 
        return x < x_jump ? eL : eR; 
    }

    /// @brief Начальная плотность
    double density(const Vector3d &r) const {
        return r.x() < x_jump ? rL : rR + epsilon * std::sin(5.0 * r.x());
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &r) const {
        return {r.x() < x_jump ? uL : uR, 0.0, 0.0};
    }

    /// @brief Начальное давление
    double pressure(const Vector3d &r) const {
        return r.x() < x_jump ? pL : pR;
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const {
        return eos.energy_rp(density(r), pressure(r));
    }

    ~ShuOsherTest() override = default;
};

}
}
