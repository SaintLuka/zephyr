#pragma once

#include <vector>

#include <zephyr/phys/fractions.h>
#include <zephyr/phys/matter/eos/eos.h>

namespace zephyr::phys {

/// @brief Дополнительные аргументы функций
struct MixOptions {
    bool deriv  = false; ///< Вычислять производные (зачем?)
    double rho0 = NAN;   ///< Начальное приближение плотности
    double P0   = NAN;   ///< Начальное приближение давления
    double T0   = NAN;   ///< Начальное приближение температуры

    /// @brief Начальные приближения для истинных плотностей
    const ScalarSet* rhos = nullptr;

    /// @brief Преобразование в опции одного УрС
    EosOptions operator[](int idx) const {
        return {.deriv = deriv,
                .rho0 = (rhos != nullptr ? (*rhos)[idx] : NAN),
                .P0 = P0, .T0 = T0};
    }
};

/// @brief Смесь с равновесием по давлению и температуре (PT-замыкание).
class MixturePT {
public:

    /// @brief Конструктор по умолчанию (пустой список)
    MixturePT(bool old_style = true) : m_old(old_style) { };

    /// @brief Создать из списка инициалиации
    MixturePT(const std::initializer_list<Eos::Ptr>& il)
        : m_materials(il) { }

    /// @brief Использовать старый метод поиска корней
    void use_old() { m_old = true; }

    /// @brief Использовать новый метод поиска корней
    void use_new() { m_old = false; }

    /// @brief Число компонент
    int size() const;

    /// @brief Удалить материалы
    void clear();

    /// @brief Добавить материал в список
    void append(Eos::Ref eos);

    /// @brief Добавить материал в список
    void operator+=(Eos::Ref eos);

    /// @brief Оператор доступа к конкретному материалу
    /// @param idx Индекс материала
    inline Eos& operator[](int idx) { return *m_materials[idx]; }

    /// @brief Оператор доступа к конкретному материалу
    /// @param idx Индекс материала
    inline const Eos& operator[](int idx) const { return *m_materials[idx]; }


    /// @brief Найти равновесное давление смеси
    /// @param density Смесевая плотность
    /// @param e Внутренняя энергия смеси
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опций целесообразно передавать начальные
    /// приближения для температуры и давления. Также указать {.deriv = true},
    /// если требуется получить производные.
    /// @details Решается метом итераций Ньютона по паре уравнений.
    dRdE pressure_re(double density, double e, const Fractions& beta,
                     const MixOptions& options = {}) const;

    /// @brief Найти равновесное давление смеси
    /// @param density Смесевая плотность
    /// @param temperature Равновесная температура
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции можно передать начальное приближение
    /// для давления.
    /// @details Решается метом итераций Ньютона по одному уравнению.
    dRdT pressure_rT(double density, double temperature,
                     const Fractions& beta, const MixOptions& = {}) const;

    /// @brief Найти внутреннюю энергию смеси
    /// @param density Смесевая плотность
    /// @param pressure Равновесное давление
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции целесообразно передавать начальное
    /// приближение для температуры.
    /// @details Решается методом итераций Ньютона по одному уравнению.
    /// Для смеси StiffenedGas точное решение получается за одну итерацию.
    dRdP energy_rP(double density, double pressure, const Fractions& beta,
                   const MixOptions& options = {}) const;

    /// @brief Вспомогательная функция, удобна для задания начальных условий.
    /// Кроме того, используется в моделях с учетом теплопроводности.
    dRdT energy_rT(double density, double temperature,
                   const Fractions& beta, const MixOptions &options = {}) const;

    /// @brief Найти равновесную температуру смеси
    /// @param density Смесевая плотность
    /// @param pressure Равновесное давление
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции можно передать начальное приближение
    /// для температуры.
    /// @details Решается метом итераций Ньютона по одному уравнению.
    /// Для смеси StiffenedGas точное решение получается за одну итерацию.
    double temperature_rP(double density, double pressure,
                          const Fractions& beta, const MixOptions& options = {}) const;

    /// @brief Смесевая скорость звука
    /// @param density Смесевая плотность
    /// @param energy Внутренняя энергия смеси
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опций целесообразно передавать начальные
    /// приближения для температуры и давления.
    /// @details Решается метом итераций Ньютона по паре уравнений.
    double sound_speed_re(double density, double energy, const Fractions& beta,
                          const MixOptions& options = {}) const;

    /// @brief Смесевая скорость звука
    /// @param density Смесевая плотность
    /// @param pressure Равновесное давление
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции можно передать начальное приближение
    /// для температуры.
    /// @details Решается метом итераций Ньютона по одному уравнению.
    /// Для смеси StiffenedGas точное решение получается за одну итерацию.
    double sound_speed_rP(double density, double pressure, const Fractions& beta,
                          const MixOptions& options = {}) const;

    /// @brief Удельный объем смеси
    /// @param pressure Равновесная плотность смеси
    /// @param temperature Равновесная температура смеси
    /// @param beta Массовые концентрации компонент
    /// @details Быстрая функция, т.к. удельный объем выражается в явном виде
    /// Очевидно, вычислительная сложность увеличивается, если в одном из
    /// уравнений состояния плотность считается неявно (частая ситуация).
    dPdT volume_PT(double pressure, double temperature,
                   const Fractions& beta, const MixOptions& = {}) const;

    /// @brief Смесевая энергия
    /// @param pressure Равновесная плотность смеси
    /// @param temperature Равновесная температура смеси
    /// @param beta Массовые концентрации компонент
    /// @details Быстрая функция, т.к. смесевая энергия выражается в явном виде.
    /// Очевидно, вычислительная сложность увеличиывется, если в одном из
    /// уравнений состояния энергия считается неявно (частая ситуация).
    dPdT energy_PT(double pressure, double temperature,
                   const Fractions& beta, const MixOptions& = {}) const;

    /// @brief Аппроксимация смеси двучленным уравнением состояния
    /// @param density Смесевая плотность
    /// @param pressure Равновесное давление
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции целесообразно передавать начальное
    /// приближение для температуры.
    /// @details Решается метом итераций Ньютона по одному уравнению.
    /// Для смеси StiffenedGas точное решение получается за одну итерацию.
    StiffenedGas stiffened_gas(double density, double pressure,
            const Fractions& beta, const MixOptions& options = {}) const;

    /// @brief Минимальное значение давления, при котором работают все
    /// уравнения состояния с ненулевыми концентрациями
    /// (максимальное из минимальных давлений)
    double min_pressure(const Fractions& beta) const;

    /// @brief Подгон теплоемкости Cv
    void adjust_cv(double rho_ref, double P_ref, double T_ref);

    /// @brief Подгон аддитивной постоянной T_0
    void adjust_T0(double rho_ref, double P_ref, double T_ref);


    // Предыдущие функции повторяют интерфейс обычных УрС.
    // Если требуется получить сразу несколько величин, то для оптимизации
    // следует вызывать функции ниже. В этом случае дорогостоящие Ньютоновские
    // итерации по системе уравнений, которые требуются для PT-замыкания,
    // будут выполнены единожды.

    using doublet_re = std::tuple<ScalarSet, dRdE>;
    using doublet_rT = std::tuple<ScalarSet, dRdT>;
    using doublet_rP = std::tuple<ScalarSet, dRdP>;

    using triplet_re = std::tuple<ScalarSet, dRdE, dRdE>;
    using triplet_rP = std::tuple<ScalarSet, dRdP, dRdP>;
    using triplet_rT = std::tuple<ScalarSet, dRdT, dRdT>;


    /// ...
    triplet_re get_rPT(double density, double energy,
                       const Fractions& beta, const MixOptions& options = {}) const;

    /// ...
    triplet_rP get_reT(double density, double pressure,
                       const Fractions& beta, const MixOptions& options = {}) const;

protected:
    // Начальное приближение для истинных плотностей
    // Если beta[i] > 0.0, тогда densities[i] адекватно определена.
    ScalarSet init_densities(const Fractions& beta, const MixOptions& options) const;



    // Метод Ньютона для поиска равновесного давления
    // (mix_density, temperature, mass_frac) -> (densities, pressure)
    // Не проверяет на чистоту
    doublet_rT find_rP_rT(double density, double temperature,
                          const Fractions& beta, const MixOptions& options = {}) const;

    // Оптимальное вычисление e(rho, T)
    triplet_rT find_reP_rT(double density, double temperature,
                           const Fractions& beta, const MixOptions& options = {}) const;

    // Метод Ньютона для поиска равновесной температуры
    // (mix_density, pressure, mass_frac) -> (vol_frac, temperature)
    // Не проверяет на чистоту
    doublet_rP find_rT_rP(double density, double pressure,
                          const Fractions& beta, const MixOptions& options = {}) const;

    // Оптимальное вычисление e(rho, P)
    triplet_rP find_reT_rP(double density, double pressure,
                           const Fractions& beta, const MixOptions& options = {}) const;

    // Метод Ньютона для поиска равновесной пары (P, T)
    // (mix_density, mix_energy, mass_frac) -> (vol_frac, pressure, temperature)
    // Не проверяет на чистоту
    triplet_re find_rPT(double density, double energy,
                        const Fractions& beta, const MixOptions& options = {}) const;


    // Старая схема. Итерации по температуре. Не проверяет на чистоту
    doublet_rT find_rP_rT_old(double density, double temperature,
                              const Fractions& beta, const MixOptions& options = {}) const;

    // Оптимальное вычисление e(rho, T)
    triplet_rT find_reP_rT_old(double density, double temperature,
                               const Fractions& beta, const MixOptions& options = {}) const;

    // Старая схема. Итерации по давлению. Не проверяет на чистоту
    doublet_rP find_rT_rP_old(double rho, double P,
                              const Fractions& beta, const MixOptions& options = {}) const;

    // Оптимальное вычисление e(rho, P)
    triplet_rP find_reT_rP_old(double rho, double P,
                               const Fractions& beta, const MixOptions& options = {}) const;

    // Старая схема. Итерации по давлению и температуре. Не проверяет на чистоту
    triplet_re find_rPT_old(double density, double e_mix,
                            const Fractions& beta, const MixOptions& options = {}) const;


    // Новая схема. Итерации по объемным долям. Не проверяет на чистоту
    doublet_rT find_rP_rT_new(double rho, double T,
                              const Fractions& beta, const MixOptions& options = {}) const;

    // Оптимальное вычисление e(rho, T)
    triplet_rT find_reP_rT_new(double density, double temperature,
                               const Fractions& beta, const MixOptions& options = {}) const;

    // Новая схема. Итерации по объемным долям и температуре. Не проверяет на чистоту
    doublet_rP find_rT_rP_new(double density, double pressure,
                              const Fractions& beta, const MixOptions& options = {}) const;

    // Оптимальное вычисление e(rho, P)
    triplet_rP find_reT_rP_new(double rho, double P,
                               const Fractions& beta, const MixOptions& options = {}) const;

    // Новая схема. Итерации по объемным долям и температуре. Не проверяет на чистоту
    triplet_re find_rPT_new(double density, double energy,
                            const Fractions& beta, const MixOptions& options = {}) const;


    /// @brief Скорость звука от равновесных температуры и давления
    double sound_speed_PT(double pressure, double temperature,
                          const Fractions& beta, const MixOptions& = {}) const;

    /// @brief Старый метод поиска корней?
    bool m_old = true;

    /// @brief Массив УрСов
    std::vector<Eos::Ptr> m_materials;
};

} // namespace zephyr::phys