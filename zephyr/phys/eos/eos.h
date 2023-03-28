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

    virtual dRdE pressure_re(double density, double energy) const;

    virtual double energy_rp(double density, double pressure) const;

    virtual double sound_speed_rp(double density, double pressure) const;


    /// @brief Прочие бесполезные формулы

    virtual double density_pe(double pressure, double energy) const;

    virtual double density_pt(double pressure, double temperature) const;

    virtual double pressure_rt(double density, double temperature) const;

    virtual double energy_rt(double density, double temperature) const;

    virtual double energy_pt(double pressure, double temperature) const;

    virtual double temperature_rp(double density, double pressure) const;

    virtual double temperature_re(double density, double energy) const;

    virtual double temperature_pe(double pressure, double energy) const;


    /// @brief Скорость звука от нетривиальных переменных

    virtual double sound_speed_re(double density, double pressure) const;

    virtual double sound_speed_rt(double density, double temperature) const;

    virtual double sound_speed_pe(double pressure, double energy) const;

    virtual double sound_speed_pt(double pressure, double temperature) const;


    /// @brief Тест для уравнений состояния, проверяет, что все формулы
    /// согласованы.
    static void test();

};

}
}
