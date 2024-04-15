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

    /// @brief Найти внутреннюю энергию смеси
    /// @param density Смесевая плотность
    /// @param pressure Равновесное давление
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции целесообразно передавать начальное
    /// приближение для температуры.
    /// @details Решается методом итераций Ньютона по одному уравнению.
    /// Для смеси StiffenedGas точное решение получается за одну итерацию.
    double energy_rp(double density, double pressure, const Fractions& beta,
                     const Options& options = {}) const;

    /// @brief Найти внутреннюю температуру и энергию смеси
    /// @param density Смесевая плотность
    /// @param pressure Равновесное давление
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции целесообразно передавать начальное
    /// приближение для температуры.
    /// @details Решается методом итераций Ньютона по одному уравнению.
    /// Для смеси StiffenedGas точное решение получается за одну итерацию.
    std::pair<double, double> temperature_energy_rp(double density, double pressure, const Fractions& beta,
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
    double sound_speed_rp(double density, double pressure, const Fractions& beta,
                          const Options& options = {}) const;

    /// @brief Найти равновесное давление смеси
    /// @param density Смесевая плотность
    /// @param temperature Равновесная температура
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции можно передать начальное приближение
    /// для давления.
    /// @details Решается метом итераций Ньютона по одному уравнению.
    double pressure_rt(double density, double temperature,
                       const Fractions& beta, const Options& = {}) const;

    /// @brief Найти равновесную температуру смеси
    /// @param density Смесевая плотность
    /// @param pressure Равновесное давление
    /// @param beta Массовые концентрации компонент
    /// @param options В качестве опции можно передать начальное приближение
    /// для температуры.
    /// @details Решается метом итераций Ньютона по одному уравнению.
    /// Для смеси StiffenedGas точное решение получается за одну итерацию.
    double temperature_rp(double density, double pressure,
                          const Fractions& beta, const Options& options = {}) const;

    /// @brief Удельный объем смеси
    /// @param pressure Равновесная плотность смеси
    /// @param temperature Равновесная температура смеси
    /// @param beta Массовые концентрации компонент
    /// @details Быстрая функция, т.к. удельный объем выражается в явном виде
    /// Очевидно, вычислительная сложность увеличивается, если в одном из
    /// уравнений состояния плотность считается неявно (частая ситуация).
    dPdT volume_pt(double pressure, double temperature, const Fractions& beta) const;

    /// @brief Смесевая энергия
    /// @param pressure Равновесная плотность смеси
    /// @param temperature Равновесная температура смеси
    /// @param beta Массовые концентрации компонент
    /// @details Быстрая функция, т.к. смесевая энергия выражается в явном виде.
    /// Очевидно, вычислительная сложность увеличиывется, если в одном из
    /// уравнений состояния энергия считается неявно (частая ситуация).
    dPdT energy_pt(double pressure, double temperature, const Fractions& beta) const;

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
PairPT find_PT(double density, double energy, const Fractions& beta,
               const Options& options = {}) const;

protected:

    /// @brief Скорость звука от равновесных температуры и давления
    double sound_speed_pt(double pressure, double temperature,
                          const Fractions& beta) const;

    /// @brief Массив материалов
    std::vector<Eos::Ptr> m_materials;
};

}