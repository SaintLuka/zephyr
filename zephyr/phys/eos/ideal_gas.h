#pragma once

#include <zephyr/phys/eos/eos.h>

namespace zephyr { namespace phys {

/// @brief Создание указателя на базовый класс Eos
/// @brief Уравнение состояния идеального газа
class IdealGas : public Eos {
public:
    double gamma;
    double Cv;

    /// @brief Конструктор с заданием параметров
    explicit IdealGas(double gamma = 1.4, double Cv = 1.0);

    /// @brief Конструктор для известных материалов
    explicit IdealGas(const std::string &name);

    /// @brief Создание указателя на базовый класс Eos
    template<typename... Args>
    static Eos::Ptr create(Args &&... args) {
        return std::make_shared<IdealGas>(std::forward<Args>(args)...);
    };

    /// @brief Основная формула. Давление от плотности и энергии
    /// @param options Передать {.deriv = true}, если необходимы производные
    dRdE pressure_re(double density, double energy,
                     const Options& options = {}) const final;

    /// @brief Основная формула. Энергия от плотности и давления
    double energy_rp(double density, double pressure,
                     const Options& = {}) const final;

    /// @brief Скорость звука от плотности и энергии
    double sound_speed_re(double density, double energy,
                          const Options& options = {}) const final;

    /// @details Скорость звука от плотности и давления
    double sound_speed_rp(double density, double pressure,
                          const Options& options = {}) const final;

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    double pressure_rt(double density, double temperature,
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

    /// @brief Аппроксимация уравнения состояния двучленным уравнением
    /// состояния в окрестности заданой плотности и давления.
    virtual StiffenedGas stiffened_gas(double density, double pressure,
                                       const Options& options = {}) const;

    /// @brief Минимальное значение давления, при котором УрС работает.
    virtual double min_pressure() const;

};

}
}
