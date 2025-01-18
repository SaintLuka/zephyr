#pragma once

#include <cmath>
#include <memory>

namespace zephyr::phys {

class StiffenedGas;

/// @brief Термодинамическая величина + производные по плотности и энергии
struct dRdE {
    double val = NAN;
    double dR = NAN;
    double dE = NAN;

    dRdE(double value, double dR = NAN, double dE = NAN)
        : val(value), dR(dR), dE(dE) { }

    /// @brief Приведение к double
    operator double() const { return val; };
};

/// @brief Термодинамическая величина + производные по плотности и температуре
struct dRdT {
    double val = NAN;
    double dR = NAN;
    double dT = NAN;

    dRdT(double value, double dR = NAN, double dT = NAN)
        : val(value), dR(dR), dT(dT) { }

    /// @brief Приведение к double
    operator double() const { return val; };
};

/// @brief Термодинамическая величина + производные по давлению и температуре
struct dPdT {
    double val = NAN;
    double dP = NAN;
    double dT = NAN;

    dPdT(double value, double dP = NAN, double dT = NAN)
        : val(value), dP(dP), dT(dT) { }

    /// @brief Приведение к double
    operator double() const { return val; };
};

/// @brief Термодинамическая величина + производные по плотности и давлению
struct dRdP {
    double val = NAN;
    double dR = NAN;
    double dP = NAN;

    dRdP(double value, double dR = NAN, double dP = NAN)
            : val(value), dR(dR), dP(dP) { }

    /// @brief Приведение к double
    operator double() const { return val; };
};

/// @brief Дополнительные аргументы функций УрС
struct EosOptions {
    bool deriv  = false; ///< Вычислять производные
    double rho0 = NAN;   ///< Начальное приближение плотности
    double P0   = NAN;   ///< Начальное приближение давления
    double T0   = NAN;   ///< Начальное приближение температуры
};

/// @brief Абстрактный класс уравнения состояния
class Eos {
public:
    /// @brief Умный указатель на базовый класс
    using Ptr = std::shared_ptr<Eos>;
    using Ref = const std::shared_ptr<Eos>&;

    /// @brief Конструктор по умолчанию
    Eos() = default;

    /// @brief Приведение к конкретному типу уравнения состояния
    /// @throw segfault при невозможности приведения
    /// @code
    ///     Eos::Ptr eos;
    ///     IdealGas& gas = eos->cast<IdealGas>();
    /// @endcode
    template <class T>
    typename std::enable_if<std::is_base_of<Eos, T>::value, T&>::type
    cast() { return *(dynamic_cast<T*>(this)); };

    /// @brief Основная формула, производные требуются для некоторых численных
    /// методов, также через первые производные можно определить скорость
    /// звука и аппроксимацию двучленным уравнением состояния.
    /// @param options Передать {.deriv = true}, если необходимы производные
    virtual dRdE pressure_re(double density, double energy,
                             const EosOptions &options = {}) const;

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    /// Кроме того, используется в моделях с учетом теплопроводности.
    virtual dRdT pressure_rT(double density, double temperature,
                             const EosOptions &options = {}) const;

    /// @brief Формула необходима, если в качестве примитивной переменной
    /// используется давление. Тогда формула позволяет вычислить энергию.
    virtual dRdP energy_rP(double density, double pressure,
                           const EosOptions &options = {}) const;

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    /// Кроме того, используется в моделях с учетом теплопроводности.
    virtual dRdT energy_rT(double density, double temperature,
                           const EosOptions &options = {}) const;

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    virtual double temperature_rP(double density, double pressure,
                                  const EosOptions &options = {}) const;

    /// @brief Скорость звука от плотности и энергии
    virtual double sound_speed_re(double density, double energy,
                                  const EosOptions &options = {}) const;

    /// @details Скорость звука от плотности и давления. При известном значении
    /// энергии целесообразно использовать функцию Eos::sound_speed_re.
    virtual double sound_speed_rP(double density, double pressure,
                                  const EosOptions &options = {}) const;

    /// @brief Удельный объем по давлению и температуре. Функция используется
    /// в формулах для PT-замыкания.
    /// @param options Для вычисления производных указать {.deriv = true},
    /// также некоторые УрС вычисляют плотность неявно, поэтому целесообразно
    /// передавать начальное приближение для плотности {.rho0 = }
    virtual dPdT volume_PT(double pressure, double temperature,
                           const EosOptions &options = {}) const;

    /// @brief Внутрення энергия по давлению и температуре. Функция
    /// используется в формулах для PT-замыкания.
    /// @param options Для вычисления производных указать {.deriv = true},
    /// также некоторые УрС внутри функции вычисляют плотность неявно, поэтому
    /// целесообразно передавать начальное приближение для плотности {.rho0 = }
    virtual dPdT energy_PT(double pressure, double temperature,
                           const EosOptions &options = {}) const;

    /// @brief Аппроксимация уравнения состояния двучленным уравнением
    /// состояния в окрестности заданой плотности и давления.
    /// В частности, функция требуется для построения численных методов
    /// на основе точного решения задачи Римана о распаде разрыва.
    virtual StiffenedGas stiffened_gas(double density, double pressure,
                                       const EosOptions &options = {}) const;

    /// @brief Референсная плотность при нормальных условиях.
    /// Точное значение не нужно, функция используется в качестве начального
    /// приближения в PT-замыкании, позволяет оценить объемные доли.
    /// В будущем можно использовать для incompressible fluids.
    double density() const { return ref_density; }

    /// @brief Минимальное давление, при котором УрС выдает приемлемые
    /// значения. При более низких значениях образуется вакуум (плотность
    /// принимает отрицательные значения, квадрат скорости звука становится
    /// отрицательным). Может принимать отрицательные значения.
    virtual double min_pressure() const;

    /// @brief Подгон теплоемкости Cv
    /// @param rho_ref, P_ref, T_ref Референсные значения плотности, давления
    /// и температуры.
    virtual void adjust_cv(double rho_ref, double P_ref, double T_ref);

    /// @brief Подгон аддитивной постоянной T_0
    /// @param rho_ref, P_ref, T_ref Референсные значения плотности, давления
    /// и температуры.
    virtual void adjust_T0(double rho_ref, double P_ref, double T_ref);

protected:
    /// @brief Референсная плотность при нормальных условиях.
    /// Точное значение не нужно, функция используется в качестве начального
    /// приближения в PT-замыкании, позволяет оценить объемные доли.
    /// В будущем можно использовать для incompressible fluids.
    double ref_density = NAN;
};

} // namespace zephyr::phys
