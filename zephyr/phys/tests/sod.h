#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/tests/test_1D.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;

/// @class Классический одномерный тест Сода
class SodTest : public Test1D {
public:
    IdealGas::Ptr eos;   ///< Используемый УрС

    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия

    /// @brief Конструктор
    SodTest() {
        eos = IdealGas::create(1.4);

        // @formatter:off
        rL = 1.0; rR = 0.125;
        pL = 1.0; pR = 0.1;
        uL = 0.0; uR = 0.0;
        eL = eos->energy_rP(rL, pL);
        eR = eos->energy_rP(rR, pR);

        x_jump = 0.5;
        finish = 0.2;
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

    std::string name() const final { return "SodTest";}

    /// @brief Левая граница области
    double xmin() const final { return 0.0; }

    /// @brief Правая граница области
    double xmax() const final { return 1.0; }

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
        return {r.x() < x_jump ? uL : uR, 0.0, 0.0};
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
        return eos;
    }
};

} // namespace zephyr::phys
