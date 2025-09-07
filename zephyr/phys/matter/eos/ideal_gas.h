#pragma once

#include <zephyr/phys/matter/eos/eos.h>

namespace zephyr::phys {

/// @brief Уравнение состояния идеального газа
class IdealGas : public Eos {
public:
    double gamma; ///< Показатель адиабаты
    double Cv;    ///< Теплоемкость при постоянном объеме

    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<IdealGas>;
    using Ref = const std::shared_ptr<IdealGas>&;

    /// @brief Конструктор с заданием параметров
    explicit IdealGas(double gamma = 1.4, double Cv = 1.0);

    /// @brief Конструктор для известных материалов
    explicit IdealGas(const std::string &name);

    /// @brief Создание указателя на базовый класс Eos
    template<typename... Args>
    static IdealGas::Ptr create(Args &&... args) {
        return std::make_shared<IdealGas>(std::forward<Args>(args)...);
    };

    /// @brief Основная формула. Давление от плотности и энергии
    /// @param options Передать {.deriv = true}, если необходимы производные
    dRdE pressure_re(double density, double energy,
                     const EosOptions& options = {}) const final;

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    dRdT pressure_rT(double density, double temperature,
                     const EosOptions& options = {}) const final;

    /// @brief Основная формула. Энергия от плотности и давления
    dRdP energy_rP(double density, double pressure,
                   const EosOptions& = {}) const final;

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    dRdT energy_rT(double density, double temperature,
                   const EosOptions& options = {}) const final;

    /// @brief Вспомогательная функция, удобна для задания начальных условий
    double temperature_rP(double density, double pressure,
                          const EosOptions& options = {}) const final;

    /// @brief Скорость звука от плотности и энергии
    double sound_speed_re(double density, double energy,
                          const EosOptions& options = {}) const final;

    /// @details Скорость звука от плотности и давления
    double sound_speed_rP(double density, double pressure,
                          const EosOptions& options = {}) const final;

    /// @brief Удельный объем по давлению и температуре. Функция используется
    /// в формулах для PT-замыкания.
    /// @param options Передать {.deriv = true}, если необходимы производные
    dPdT volume_PT(double pressure, double temperature,
                   const EosOptions& options = {}) const final;

    /// @brief Внутренняя энергия по давлению и температуре. Функция
    /// используется в формулах для PT-замыкания.
    /// @param options Передать {.deriv = true}, если необходимы производные
    dPdT energy_PT(double pressure, double temperature,
                   const EosOptions& options = {}) const final;

    /// @brief Аппроксимация уравнения состояния двучленным уравнением
    /// состояния в окрестности заданной плотности и давления.
    StiffenedGas stiffened_gas(double density, double pressure,
                               const EosOptions& options = {}) const final;

    /// @brief Минимальное значение давления, при котором УрС работает
    double min_pressure() const final;

    /// @brief Подгон теплоемкости Cv
    void adjust_cv(double density, double pressure, double temperature) final;

};

} // namespace zephyr::phys
