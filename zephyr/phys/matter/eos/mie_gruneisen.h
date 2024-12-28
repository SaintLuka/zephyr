#pragma once

#include <zephyr/phys/matter/eos/eos.h>

namespace zephyr::phys {

/// @brief Непонятное уравнение состояния, какая-то смесь уравнений
/// состояния Ми-Грюнайзена и Мурнагана, хорошо было бы найти исходник.
class MieGruneisen : public Eos {
public:
    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<MieGruneisen>;

    /// @brief Конструктор для известных материалов
    explicit MieGruneisen(const std::string &name);

    /// @brief Создание указателя на базовый класс Eos
    template<typename... Args>
    static MieGruneisen::Ptr create(Args &&... args) {
        return std::make_shared<MieGruneisen>(std::forward<Args>(args)...);
    };

    /// @brief Основная формула. Давление от плотности и энергии
    /// @param options Передать {.deriv = true}, если необходимы производные
    dRdE pressure_re(double density, double energy,
                     const Options& options = {}) const final;

    /// @brief Основная формула. Энергия от плотности и давления
    double energy_rP(double density, double pressure,
                     const Options& options = {}) const final;

    /// @brief Скорость звука от плотности и энергии
    double sound_speed_re(double density, double energy,
                          const Options& options = {}) const final;

    /// @details Скорость звука от плотности и давления
    double sound_speed_rP(double density, double pressure,
                          const Options& options = {}) const final;

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    dRdT pressure_rT(double density, double temperature,
                     const Options& options = {}) const final;

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    /// Кроме того, используется в моделях с учетом теплопроводности.
    dRdT energy_rT(double density, double temperature,
                   const Options &options = {}) const final;

    /// @brief Вспомогательная функция, удобна для задания начальных условий
    double temperature_rP(double density, double pressure,
                          const Options& options = {}) const final;

    /// @brief Удельный объем по давлению и температуре. Функция используется
    /// в формулах для PT-замыкания.
    /// @param options Для вычисления производных указать {.deriv = true},
    /// Данный УрС вычисляет плотность неявно, поэтому целесообразно
    /// передавать начальное приближение для плотности {.rho0 = }
    dPdT volume_PT(double pressure, double temperature,
                   const Options& options = {}) const final;

    /// @brief Внутренняя энергия по давлению и температуре. Функция
    /// используется в формулах для PT-замыкания.
    /// @param options Для вычисления производных указать {.deriv = true},
    /// Для данного уравнения состояния в функции выполняется неявная процедура
    /// для нахождения плотнсти, поэтому целесообразно передавать начальное
    /// приближение для плотности {.rho0 = }
    dPdT energy_PT(double pressure, double temperature,
                   const Options& options = {}) const final;

    /// @brief Аппроксимация уравнения состояния двучленным уравнением
    /// состояния в окрестности заданой плотности и давления.
    /// @param options
    virtual StiffenedGas stiffened_gas(double density, double pressure,
                                       const Options& options = {}) const;

    /// @brief Минимальное значение давления, при котором УрС работает.
    virtual double min_pressure() const;

public:

    double cold_pressure(double density) const;

    double cold_energy(double density) const;

    /// @brief Задать параметры материала из таблицы
    void table_params(std::string name);


    double v0; ///< Референсный удельный объем
    double B;  ///< Референсный модуль упругости
    double n;  ///< Безразмерный коэффициент (степень)
    double Gr; ///< Коэффициент Грюнайзена

    double T0; ///< Референсная температура
    double Cv; ///< Теплоемкость при постоянном давлении
};

} // namespace zephyr::phys
