#pragma once

#include <zephyr/phys/eos/eos.h>

namespace zephyr { namespace phys {

class IdealGas : public Eos {
public:
    double gamma;
    double Cv;

    /// @brief Конструктор с заданием параметров
    explicit IdealGas(double gamma = 1.4, double Cv = 1.0);

    /// @brief Конструктор для известных материалов
    explicit IdealGas(const std::string &name);

    template<typename... Args>
    static Eos::Ptr create(Args &&... args) {
        return std::make_shared<IdealGas>(std::forward<Args>(args)...);
    };


    /// @brief Основные формулы для уравнений состояния, необходимые для
    /// решения уравнений Эйлера. Только связь внутренней энергии, плотности
    /// и давления

    dRdE pressure_re(double density, double energy, const Options& = {}) const final;

    double energy_rp(double density, double pressure, const Options& = {}) const final;

    double sound_speed_re(double density, double energy, const Options& = {}) const final;

    double sound_speed_rp(double density, double pressure, const Options& = {}) const final;


    /// @brief Удобно для задания начальных условий

    double pressure_rt(double density, double temperature, const Options& = {}) const final;

    double temperature_rp(double density, double pressure, const Options& = {}) const final;


    /// @brief Следующие функции используются для PT-замыкания

    dPdT volume_pt(double pressure, double temperature, const Options& = {}) const final;

    dPdT energy_pt(double pressure, double temperature, const Options& = {}) const final;


    /// @brief Аппроксимация двучленным УрС

    virtual StiffenedGas stiffened_gas(double density, double pressure, const Options& = {}) const;

};

}
}
