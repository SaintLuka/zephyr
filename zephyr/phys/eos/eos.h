#pragma once

#include <memory>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/eos/types.h>

namespace zephyr { namespace phys {

class StiffenedGas;

/// @brief Абстрактный класс уравнения состояния
class Eos {
public:
    /// @brief Умный указатель на базовый класс
    using Ptr = std::shared_ptr<Eos>;

    /// @brief Конструктор по умолчанию
    Eos() = default;

    /// @brief Основные формулы для уравнений состояния, необходимые для
    /// решения уравнений Эйлера. Только связь внутренней энергии, плотности
    /// и давления

    virtual dRdE pressure_re(double density, double energy, const Options& = {}) const;

    virtual double energy_rp(double density, double pressure, const Options& = {}) const;

    virtual double sound_speed_re(double density, double energy, const Options& = {}) const;

    /// @details При известной энергии целесообразно использовать функцию
    /// sound_speed_re.
    virtual double sound_speed_rp(double density, double pressure, const Options& = {}) const;


    /// @brief Удобно для задания начальных условий

    virtual double pressure_rt(double density, double temperature, const Options& = {}) const;

    virtual double temperature_rp(double density, double pressure, const Options& = {}) const;


    /// @brief Следующие функции используются для PT-замыкания

    virtual dPdT volume_pt(double pressure, double temperature, const Options& = {}) const;

    virtual dPdT energy_pt(double pressure, double temperature, const Options& = {}) const;


    /// @brief Аппроксимация двучленным УрС

    virtual StiffenedGas stiffened_gas(double density, double pressure, const Options& = {}) const;


    /// @brief Тест для уравнений состояния, проверяет, что все формулы
    /// согласованы.
    static void test();

};

}
}
