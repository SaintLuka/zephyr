#pragma once

#include <zephyr/phys/matter/eos/eos.h>
#include <zephyr/phys/matter/cond/conductivity.h>
#include <zephyr/phys/matter/visc/viscosity.h>
#include <zephyr/phys/matter/plast/plasticity.h>

namespace zephyr::phys {

class StiffenedGas;

/// @brief Класс материала, содержит различные свойства, как минимум
/// обязан содержать уравнение состояния.
class Material {
public:
    /// @brief Создать материал по уравнению состояния
    template<class SomeEos>
    Material(std::shared_ptr<SomeEos> eos) {
        static_assert(std::is_base_of<Eos, SomeEos>::value,
                      "Can create material only from EoS");

        m_eos   = eos;
        m_cond  = Conductivity::create();
        m_visc  = Viscosity::create();
        m_plast = Plasticity::create();
    }

    /// @brief Тэг материала
    const std::string &tag() const;

    /// @brief Задать тэг материала
    void set_tag(const std::string &tag);


    //  Добавление (или замена) свойств материала

    /// @brief Добавить (заменить) уравнение состояния
    void add(Eos::Ref eos);

    /// @brief Добавить (заменить) модель теплопроводности
    void add(Conductivity::Ref cond);

    /// @brief Добавить (заменить) модель вязкости
    void add(Viscosity::Ref visc);

    /// @brief Добавить (заменить) модель пластичности
    void add(Plasticity::Ref plast);

    /// @brief Добавить (заменить) уравнение состояния
    void operator+=(Eos::Ref eos);

    /// @brief Добавить (заменить) модель теплопроводности
    void operator+=(Conductivity::Ref cond);

    /// @brief Добавить (заменить) модель вязкости
    void operator+=(Viscosity::Ref visc);

    /// @brief Добавить (заменить) модель пластичности
    void operator+=(Plasticity::Ref plast);


    // Вытащить различные свойства из материала (обычно не нужно)

    /// @brief Получить уравнение состояния
    Eos::Ref eos() const { return m_eos; }

    /// @brief Получить модель теплопроводности
    Conductivity::Ref cond() const { return m_cond; }

    /// @brief Получить модель вязкости
    Viscosity::Ref visc() const { return m_visc; }

    /// @brief Получить модель пластичности
    Plasticity::Ref plast() const { return m_plast; }


    // Функции, связанные с уравнением состояния

    /// @brief Референсная плотность, плотность для Incompressible
    double density() const;

    /// @brief Давление от плотности и внутренней энергии
    dRdE pressure_re(double density, double energy, const EosOptions &options = {}) const;

    /// @brief Давление от плотности и температуры
    dRdT pressure_rT(double density, double temperature, const EosOptions &options = {}) const;

    /// @brief Внутренняя энергия от плотности и температуры
    dRdT energy_rT(double density, double temperature, const EosOptions &options = {}) const;

    /// @brief Скорость звука от плотности и внутренней энергии
    double sound_speed_re(double density, double energy, const EosOptions &options = {}) const;

    /// @details Скорость звука от плотности и давления
    double sound_speed_rP(double density, double pressure, const EosOptions &options = {}) const;

    /// @brief Внутренняя энергия от плотности и давления
    double energy_rP(double density, double pressure, const EosOptions &options = {}) const;

    /// @brief Температура от плотности и давления
    double temperature_rP(double density, double pressure, const EosOptions &options = {}) const;

    /// @brief Удельный объем по давлению и температуре.
    dPdT volume_PT(double pressure, double temperature, const EosOptions &options = {}) const;

    /// @brief Внутрення энергия по давлению и температуре.
    dPdT energy_PT(double pressure, double temperature, const EosOptions &options = {}) const;

    /// @brief Аппроксимация двучленным уравнением состояния в окрестности заданной плотности и давления.
    StiffenedGas stiffened_gas(double density, double pressure, const EosOptions &options = {}) const;

    /// @brief Минимальное давление, при котором УрС выдает приемлемые значения.
    double min_pressure() const;

    /// @brief Подгон теплоемкости Cv
    void adjust_cv(double rho_ref, double P_ref, double T_ref);

    /// @brief Подгон аддитивной постоянной T_0
    void adjust_T0(double rho_ref, double P_ref, double T_ref);


    // Теплопроводность

    /// @brief Коэффициент теплопроводности
    double kappa(double temperature = 0.0) const;


    // Вязкость

    /// @brief Кинематическая вязкость
    double kinematic_visc(double temperature = 0.0) const;

    /// @brief Свдиговая/динамическая вязкость
    double shear_visc(double temperature = 0.0) const;

    /// @brief Объемная/вторая вязкость
    double volume_visc(double temperature = 0.0) const;


    // Упруго-пластичность

    /// @brief Модуль сдвига
    double shear() const;

    /// @brief Модуль Юнга
    double young() const;

    /// @brief Коэффициент Пуассона
    double poisson() const;

    /// @brief Предел текучести
    double yield() const;


private:
    /// @brief Тэг материала, может быть не связан с названием материала,
    /// а представлять роль в расчете, к примеру.
    std::string m_tag;

    Eos::Ptr          m_eos;    ///< Уравнение состояния
    Conductivity::Ptr m_cond;   ///< Теплопроводность
    Viscosity::Ptr    m_visc;   ///< Вязкость
    Plasticity::Ptr   m_plast;  ///< Пластичность
};

} // namespace zephyr::phys