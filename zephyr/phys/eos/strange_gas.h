#pragma once

#include <zephyr/phys/eos/eos.h>

namespace zephyr::phys {

/// @brief Двучленное уравнение состояния
class StrangeGas : public Eos {
public:
    double gamma;
    double Cv;
    double P0;
    double eps_0;
    double T0;
    double Pe;

    /// @brief Конструктор с заданием параметров
    explicit StrangeGas(double gamma, double p_inf = 0.0,
                          double eps_0 = 0.0, double Cv = 1.0, double T0 = 0.0);

    /// @brief Конструктор для известных материалов
    explicit StrangeGas(const std::string &name);

    /// @brief Создание указателя на базовый класс Eos
    template<typename... Args>
    static Eos::Ptr create(Args &&... args) {
        return std::make_shared<StrangeGas>(std::forward<Args>(args)...);
    };

    /// @brief Основная формула. Давление от плотности и энергии
    /// @param options Передать {.deriv = true}, если необходимы производные
    dRdE pressure_re(double density, double energy,
                     const Options& options = {}) const final;

    /// @brief Основная формула. Энергия от плотности и давления
    double energy_rp(double density, double pressure,
                     const Options& options = {}) const final;
/*
    /// @brief Скорость звука от плотности и энергии
    double sound_speed_re(double density, double energy,
                          const Options& options = {}) const final;

    /// @details Скорость звука от плотности и давления
    double sound_speed_rp(double density, double pressure,
                          const Options& options = {}) const final;
*/
    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    dRdT pressure_rt(double density, double temperature,
                     const Options& options = {}) const final;

    /// @brief Вспомогательная функция, удобна для задания начальных условий
    double temperature_rp(double density, double pressure,
                          const Options& options = {}) const final;

    /// @brief Удельный объем по давлению и температуре. Функция используется
    /// в формулах для PT-замыкания.
    /// @param options Передать {.deriv = true}, если необходимы производные
    dPdT volume_pt(double pressure, double temperature,
                   const Options& options = {}) const final;

    /// @brief Внутренняя энергия по давлению и температуре. Функция
    /// используется в формулах для PT-замыкания.
    /// @param options Передать {.deriv = true}, если необходимы производные
    dPdT energy_pt(double pressure, double temperature,
                   const Options& options = {}) const final;
/*
    /// @brief Аппроксимация уравнения состояния двучленным уравнением
    /// состояния в окрестности заданой плотности и давления.
    virtual StiffenedGas stiffened_gas(double density, double pressure,
                                       const Options& options = {}) const;
*/
    /// @brief Минимальное значение давления, при котором УрС работает.
    virtual double min_pressure() const;
    
protected:

    double P_fix(double P) const;
    
    double dP_fix(double P) const;

    double P_inv(double P_fixed) const;

    double dP_inv(double P_fixed) const;

};

}
