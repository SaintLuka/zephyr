#pragma once

#include <zephyr/phys/literals.h>
#include <zephyr/phys/eos/types.h>

namespace zephyr { namespace phys {

/// @brief Абстрактный класс уравнения состояния
class Eos {
public:

    Eos() = default;

    /// @brief Основные формулы для уравнений состояния, необходимые для
    /// решения уравнений Эйлера. Только связь внутренней энергии, плотности
    /// и давления

    virtual dRdE pressure_re(double density, double energy, const Options& = {}) const;

    virtual dRdP energy_rp(double density, double pressure, const Options& = {}) const;

    virtual double sound_speed_re(double density, double energy, const Options& = {}) const;

    virtual double sound_speed_rp(double density, double pressure, const Options& = {}) const;


    /// @brief Удобно для задания начальных условий

    virtual double pressure_rt(double density, double temperature, const Options& = {}) const;


    /// @brief Следующие функции используются для PT-замыкания

    virtual dPdT density_pt(double pressure, double temperature, const Options& = {}) const;

    virtual dPdT energy_pt(double pressure, double temperature, const Options& = {}) const;

    virtual double temperature_rp(double density, double pressure, const Options& = {}) const;


    /// @brief Аппроксимация двучленным УрС

    virtual double stiff_gamma(double density, double pressure, const Options& = {}) const;

    virtual double stiff_p0(double density, double pressure, const Options& = {}) const;

    virtual double stiff_e0(double density, double pressure, const Options& = {}) const;


    /// @brief Тест для уравнений состояния, проверяет, что все формулы
    /// согласованы.
    static void test();

};

}
}
