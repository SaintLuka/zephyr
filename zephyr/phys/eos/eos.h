#pragma once

#include <memory>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/eos/types.h>

namespace zephyr::phys {

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

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    /// Кроме того, используется в моделях с учетом теплопроводности.
    virtual dRdT pressure_rT(double density, double temperature,
                             const Options &options = {}) const;

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    /// Кроме того, используется в моделях с учетом теплопроводности.
    virtual dRdT energy_rT(double density, double temperature,
                             const Options &options = {}) const;

    /// @brief Скорость звука от плотности и энергии
    virtual double sound_speed_re(double density, double energy,
                                  const Options &options = {}) const;

    /// @details Скорость звука от плотности и давления. При известном значении
    /// энергии целесообразно использовать функцию Eos::sound_speed_re.
    virtual double sound_speed_rP(double density, double pressure,
                                  const Options &options = {}) const;

    /// @brief Формула необходима, если в качестве примитивной переменной
    /// используется давление. Тогда формула позволяет вычислить энергию.
    virtual double energy_rP(double density, double pressure,
                             const Options &options = {}) const;

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    virtual double temperature_rP(double density, double pressure,
                                  const Options &options = {}) const;

    /// @brief Удельный объем по давлению и температуре. Функция используется
    /// в формулах для PT-замыкания.
    /// @param options Для вычисления производных указать {.deriv = true},
    /// также некоторые УрС вычисляют плотность неявно, поэтому целесообразно
    /// передавать начальное приближение для плотности {.rho0 = }
    virtual dPdT volume_PT(double pressure, double temperature,
                           const Options &options = {}) const;

    /// @brief Внутрення энергия по давлению и температуре. Функция
    /// используется в формулах для PT-замыкания.
    /// @param options Для вычисления производных указать {.deriv = true},
    /// также некоторые УрС внутри функции вычисляют плотность неявно, поэтому
    /// целесообразно передавать начальное приближение для плотности {.rho0 = }
    virtual dPdT energy_PT(double pressure, double temperature,
                           const Options &options = {}) const;

    /// @brief Аппроксимация уравнения состояния двучленным уравнением
    /// состояния в окрестности заданой плотности и давления.
    /// В частности, функция требуется для построения численных методов
    /// на основе точного решения задачи Римана о распаде разрыва.
    virtual StiffenedGas stiffened_gas(double density, double pressure,
                                       const Options &options = {}) const;

    /// @brief Минимальное давление, при котором УрС выдает приемлемые
    /// значения. При более низких значениях образуется вакуум (плотность
    /// принимает отрицательные значения, квадрат скорости звука становится
    /// отрицательным). Может принимать отрицательные значения.
    virtual double min_pressure() const;

    /// @brief Референсная плотность, значение используется в качестве
    /// начального приближения в PT-замыкании.
    double ref_density() const;

    /// @brief Подгон теплоемкости Cv
    /// @param rho_ref, P_ref, T_ref Референсные значения плотности, давления
    /// и температуры.
    virtual void adjust_cv(double rho_ref, double P_ref, double T_ref);

protected:
    /// @brief Референсная плотность при нормальных условиях.
    /// Точное значение не нужно, значение используется в качестве
    /// начального приближения в PT-замыкании.
    double rho_0 = 1.0;
};

} // namespace zephyr::phys
