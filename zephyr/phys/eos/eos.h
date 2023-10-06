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

    /// @brief Основная формула, производные требуются для некоторых численных
    /// методов, также через первые производные можно определить скорость
    /// звука и аппроксимацию двучленным уравнением состояния.
    /// @param options Передать {.deriv = true}, если необходимы производные
    virtual dRdE pressure_re(double density, double energy,
                             const Options &options = {}) const;

    /// @brief Формула необходима, если в качестве примитивной переменной
    /// используется давление. Тогда формула позволяет вычислить энергию.
    virtual double energy_rp(double density, double pressure,
                             const Options &options = {}) const;

    /// @brief Скорость звука от плотности и энергии
    virtual double sound_speed_re(double density, double energy,
                                  const Options &options = {}) const;

    /// @details Скорость звука от плотности и давления. При известном значении
    /// энергии целесообразно использовать функцию Eos::sound_speed_re.
    virtual double sound_speed_rp(double density, double pressure,
                                  const Options &options = {}) const;

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    /// Кроме того, используется в моделях с учетом теплопроводности.
    virtual double pressure_rt(double density, double temperature,
                               const Options &options = {}) const;

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    virtual double temperature_rp(double density, double pressure,
                                  const Options &options = {}) const;

    /// @brief Удельный объем по давлению и температуре. Функция используется
    /// в формулах для PT-замыкания.
    /// @param options Для вычисления производных указать {.deriv = true},
    /// также некоторые УрС вычисляют плотность неявно, поэтому целесообразно
    /// передавать начальное приближение для плотности {.rho0 = }
    virtual dPdT volume_pt(double pressure, double temperature,
                           const Options &options = {}) const;

    /// @brief Внутрення энергия по давлению и температуре. Функция
    /// используется в формулах для PT-замыкания.
    /// @param options Для вычисления производных указать {.deriv = true},
    /// также некоторые УрС внутри функции вычисляют плотность неявно, поэтому
    /// целесообразно передавать начальное приближение для плотности {.rho0 = }
    virtual dPdT energy_pt(double pressure, double temperature,
                           const Options &options = {}) const;

    /// @brief Аппроксимация уравнения состояния двучленным уравнением
    /// состояния в окрестности заданой плотности и давления.
    /// В частности, функция требуется для построения численных методов
    /// на основе точного решения задачи Римана о распаде разрыва.
    virtual StiffenedGas stiffened_gas(double density, double pressure,
                                       const Options &options = {}) const;

};

}
}
