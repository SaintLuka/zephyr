#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/eos/stiffened_gas.h>
#include <zephyr/phys/tests/test_1D.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;

/// @brief Тест с разреженной водой.
/// Отрицательное начальное давление воды, контакт с газом.
class RarefiedWater : public Test1D {
public:
    StiffenedGas::Ptr water;  ///< УрС воды
    IdealGas::Ptr     air;    ///< УрС газа

    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия

    /// @brief Конструктор
    RarefiedWater() {
        water = StiffenedGas::create("Water");
        air   = IdealGas::create("Air");

        double T = 300.0;

        rL = 0.9_g_cm3;
        uL = 0.0;
        pL = water->pressure_rT(rL, T);
        eL = water->energy_rP(rL, pL);

        rR = 1.16_kg_m3;
        uR = 0.0;
        pR = air->pressure_rT(rR, T);
        eR = air->energy_rP(rR, pR);

        x_jump = 0.9_cm;
        finish = 5.0_us;
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

    std::string name() const final { return "Rarefied water";}

    /// @brief Левая граница области
    double xmin() const final { return 0.0; }

    /// @brief Правая граница области
    double xmax() const final { return 1.0_cm; }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    ///@brief Получить положение разрыва
    double get_x_jump() const final { return x_jump; }


    /// @brief Начальная плотность
    double density(const Vector3d& r) const final {
        return r.x() < x_jump ? rL : rR;
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d& r) const final {
        return { r.x() < x_jump ? uL : uR, 0, 0};
    }

    /// @brief Начальное давление
    double pressure(const Vector3d& r) const final {
        return r.x() < x_jump ? pL : pR;
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d& r) const final {
        return r.x() < x_jump ? eL : eR;
    }

    /// @brief Уравнение состояния
    Eos::Ptr get_eos(const Vector3d& r) const final {
        return r.x() < x_jump ? (Eos::Ptr) water : (Eos::Ptr) air;
    }
};

} // namespace zephyr::phys
