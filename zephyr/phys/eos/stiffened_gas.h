#pragma once

#include <zephyr/phys/eos/eos.h>

namespace zephyr { namespace phys {

class StiffenedGas : public Eos {
public:
    double gamma;
    double Cv;
    double p_inf;
    double eps_0;

    /// @brief Конструктор для известных материалов
    StiffenedGas(const std::string &name);


    /// @brief Основные формулы для уравнений состояния, необходимые для
    /// решения уравнений Эйлера. Только связь внутренней энергии, плотности
    /// и давления

    dRdE pressure_re(double density, double energy) const final;

    double energy_rp(double density, double pressure) const final;

    double sound_speed_rp(double density, double pressure) const final;


    /// @brief Прочие бесполезные формулы

    double density_pe(double pressure, double energy) const final;

    double density_pt(double pressure, double temperature) const final;

    double pressure_rt(double density, double temperature) const final;

    double energy_rt(double density, double temperature) const final;

    double energy_pt(double pressure, double temperature) const final;

    double temperature_rp(double density, double pressure) const final;

    double temperature_re(double density, double energy) const final;

    double temperature_pe(double pressure, double energy) const final;


    /// @brief Скорость звука от нетривиальных переменных

    double sound_speed_re(double density, double pressure) const final;

    double sound_speed_rt(double density, double temperature) const final;

    double sound_speed_pe(double pressure, double energy) const final;

    double sound_speed_pt(double pressure, double temperature) const final;

};

}
}
