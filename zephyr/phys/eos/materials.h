#pragma once

#include <vector>

#include <zephyr/phys/fractions.h>
#include <zephyr/phys/eos/eos.h>
#include <zephyr/phys/eos/stiffened_gas.h>


namespace zephyr::phys {

/// @brief Набор материалов.
/// Можно использовать в качестве смесевого УрС (PT-замыкание).
class Materials {
public:

    /// @brief Конструктор по умолчанию (пустой список)
    Materials() = default;

    /// @brief Число компонент
    int size() const;

    /// @brief Удалить материалы
    void clear();

    /// @brief Добавить материал в список
    void append(Eos::Ptr eos);

    /// @brief Добавить материал в список
    void operator+=(Eos::Ptr eos);

    /// @brief Оператор доступа к конкретному материалу
    /// @param idx Индекс материала
    Eos& operator[](int idx);

    /// @brief Оператор доступа к конкретному материалу
    /// @param idx Индекс материала
    const Eos& operator[](int idx) const;

    /// @brief Найти равновесное давление смеси
    /// @param density Смесевая плотность
    /// @param energy Внутренняя энергия смеси
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опций целесообразно передавать начальные
    /// приближения для температуры и давления. Также указать {.deriv = true},
    /// если требуется получить производные.
    /// @details Решается метом итераций Ньютона по паре уравнений.
    dRdE pressure_re(double density, double energy, const Fractions& beta,
                     const Options& options = {}) const;

    dRdE pressure_re2(double density, double energy, const Fractions& beta,
                     const Options& options = {}) const;

    /// @brief Найти внутреннюю энергию смеси
    /// @param density Смесевая плотность
    /// @param pressure Равновесное давление
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции целесообразно передавать начальное
    /// приближение для температуры.
    /// @details Решается методом итераций Ньютона по одному уравнению.
    /// Для смеси StiffenedGas точное решение получается за одну итерацию.
    double energy_rP(double density, double pressure, const Fractions& beta,
                     const Options& options = {}) const;

    /// @brief Смесевая скорость звука
    /// @param density Смесевая плотность
    /// @param energy Внутренняя энергия смеси
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опций целесообразно передавать начальные
    /// приближения для температуры и давления.
    /// @details Решается метом итераций Ньютона по паре уравнений.
    double sound_speed_re(double density, double energy, const Fractions& beta,
                          const Options& options = {}) const;

    /// @brief Смесевая скорость звука
    /// @param density Смесевая плотность
    /// @param pressure Равновесное давление
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции можно передать начальное приближение
    /// для температуры.
    /// @details Решается метом итераций Ньютона по одному уравнению.
    /// Для смеси StiffenedGas точное решение получается за одну итерацию.
    double sound_speed_rP(double density, double pressure, const Fractions& beta,
                          const Options& options = {}) const;

    /// @brief Найти равновесное давление смеси
    /// @param density Смесевая плотность
    /// @param temperature Равновесная температура
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции можно передать начальное приближение
    /// для давления.
    /// @details Решается метом итераций Ньютона по одному уравнению.
    double pressure_rT(double density, double temperature,
                       const Fractions& beta, const Options& = {}) const;

    /// @brief Найти равновесную температуру смеси
    /// @param density Смесевая плотность
    /// @param pressure Равновесное давление
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции можно передать начальное приближение
    /// для температуры.
    /// @details Решается метом итераций Ньютона по одному уравнению.
    /// Для смеси StiffenedGas точное решение получается за одну итерацию.
    double temperature_rP(double density, double pressure,
                          const Fractions& beta, const Options& options = {}) const;

    /// @brief Удельный объем смеси
    /// @param pressure Равновесная плотность смеси
    /// @param temperature Равновесная температура смеси
    /// @param beta Массовые концентрации компонент
    /// @details Быстрая функция, т.к. удельный объем выражается в явном виде
    /// Очевидно, вычислительная сложность увеличивается, если в одном из
    /// уравнений состояния плотность считается неявно (частая ситуация).
    dPdT volume_PT(double pressure, double temperature,
                   const Fractions& beta, const Options& = {}) const;

    /// @brief Смесевая энергия
    /// @param pressure Равновесная плотность смеси
    /// @param temperature Равновесная температура смеси
    /// @param beta Массовые концентрации компонент
    /// @details Быстрая функция, т.к. смесевая энергия выражается в явном виде.
    /// Очевидно, вычислительная сложность увеличиывется, если в одном из
    /// уравнений состояния энергия считается неявно (частая ситуация).
    dPdT energy_PT(double pressure, double temperature,
                   const Fractions& beta, const Options& = {}) const;

    /// @brief Аппроксимация смеси двучленным уравнением состояния
    /// @param density Смесевая плотность
    /// @param pressure Равновесное давление
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции целесообразно передавать начальное
    /// приближение для температуры.
    /// @details Решается метом итераций Ньютона по одному уравнению.
    /// Для смеси StiffenedGas точное решение получается за одну итерацию.
    StiffenedGas stiffened_gas(double density, double pressure,
            const Fractions& beta, const Options& options = {}) const;

    /// @brief Минимальное значение давления, при котором работают все
    /// уравнения состояния с ненулевыми концентрациями
    /// (максимальное из минимальных давлений)
    double min_pressure(const Fractions& beta) const;

    /// @brief Найти равновесные температуру и давление смеси
    /// @param options В качестве опций целесообразно передавать начальные
    /// приближения для температуры и давления.
    /// @details Решается метом итераций Ньютона по паре уравнений.
    PairPT find_PT(double density, double energy,
                   const Fractions& beta, const Options& options = {}) const;

    /// @brief Найти равновесные температуру и внутреннюю энергию смеси
    /// @param options В качестве опций целесообразно передавать начальные
    /// приближение для температуры.
    PairET find_eT(double density, double pressure,
                   const Fractions& beta, const Options& options = {}) const;

protected:
    /// Классическая, метод Ньютона, использует функции volume(P, T)
    double pressure_rT_ver1(double density, double temperature,
                            const Fractions& beta, const Options& = {}) const;

    // Метод Ньютона, использует функции P(rho, T)
    double pressure_rT_ver2(double density, double temperature,
                            const Fractions& beta, const Options& = {}) const;

    /// Классическая, метод Ньютона, использует функции
    ///   volume(P, T), energy(P, T)
    PairPT find_PT_ver1(double density, double energy,
                        const Fractions& beta, const Options& options = {}) const;

    PairPT find_PT_ver2(double density, double energy,
                        const Fractions& beta, const Options& options = {}) const;

    /// @brief Обновить объемные доли (alpha) на новом приближении
    /// итерационного алгоритма
    /// @param rho Плотность смеси
    /// @param P Давление (близкое к равновесию)
    /// @param T Температура (близкая к равновесию)
    /// @param beta Массовые концентрации
    /// @param alpha Объемные доли, содержат приближение и обновляютя
    void update_alpha(double rho, double P, double T,
                      const Fractions& beta, Fractions& alpha) const;

    /// @brief Скорость звука от равновесных температуры и давления
    double sound_speed_PT(double pressure, double temperature,
                          const Fractions& beta, const Options& = {}) const;

    /// @brief Массив материалов
    std::vector<Eos::Ptr> m_materials;
};

}