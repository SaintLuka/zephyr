#pragma once

#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/matter/eos/stiffened_gas.h>

namespace zephyr::math {

using phys::StiffenedGas;

/// @brief Класс с набором функций для получения точного решения 
/// задачи Римана о распаде разрыва. 
/// @details Наиболее эффективные (быстрые) функции являются статичными,
/// экземпляр класса требуется создавать в том случае, если необходимо
/// получить полное решение с зависимостью от времени.
class RiemannSolver {
public:
    using cref = const double&;

    /// @brief Минималистичная структура для хранения решения
    /// задачи Римана о распаде разрыва на грани.
    /// Структура содержит только те величины, которые необходимы для
    /// вычисления потока в одноматериальной газодинамической задаче.
    struct Solution {
        double rho; ///< Плотность на грани
        double Uf;  ///< Скорость на грани
        double Pf;  ///< Давление на грани
        double P;   ///< Давление на контактном разрыве
        double U;   ///< Скорость контактного разрыва
        bool conv;  ///< Хорошая сходимость

        /// @brief Конструктор по умолчанию
        Solution() : rho(0.0 / 0.0), Uf(0.0 / 0.0), Pf(0.0 / 0.0), conv(false) {}

        /// @brief Простейший конструктор
        Solution(cref rho, cref Uf, cref Pf, cref U, cref  P, bool conv)
            : rho(rho), Uf(Uf), Pf(Pf), U(U), P(P), conv(conv) {}
    };

    /// @brief Инициализирует NAN
    RiemannSolver();

    /// @brief Решить задачу Римана о распаде разрыва
    /// @param rhoL, uL, pL Плотность, скорость, давление слева
    /// @param rhoR, uR, pR Плотность, скорость, давление справа
    /// @param gL, p0L, e0L Параметра материала слева (двучленный УрС)
    /// @param gR, p0R, e0R Параметры материала справа (двучленный УрС)
    /// @return Решение на грани
    static Solution solve(
            cref rL, cref uL, cref pL, cref gL, cref p0L,
            cref rR, cref uR, cref pR, cref gR, cref p0R);

    /// @brief Решить задачу Римана о распаде разрыва
    /// @param zL Вектор состояния слева
    /// @param zR Вектор состояния справа
    /// @param eos Уравнение состояния
    /// @return Решение на грани
    static Solution solve(const smf::PState &zL,
            const smf::PState &zR, const StiffenedGas &eos);

    /// @brief Решить задачу Римана о распаде разрыва
    /// @param zL Вектор состояния слева
    /// @param zR Вектор состояния справа
    /// @param eosL Уравнение состояния слева
    /// @param eosR Уравнение состояния справа
    /// @return Решение на грани
    static Solution solve(const smf::PState &zL, const smf::PState &zR,
            const StiffenedGas &eosL, const StiffenedGas &eosR);

    /// @brief Одноматериальный конструктор
    /// @param zL Вектор состояния слева
    /// @param zR Вектор состояния справа
    /// @param eos Уравнение состояния
    /// @param x_jump Положение разрыва
    RiemannSolver(const smf::PState &zL, const smf::PState &zR,
                  const StiffenedGas &eos, double x_jump = 0.0);

    /// @brief Двухматериальный конструктор
    /// @param zL Вектор состояния слева
    /// @param zR Вектор состояния справа
    /// @param eosL Уравнение состояния слева
    /// @param eosR Уравнение состояния справа
    /// @param x_jump Положение разрыва
    RiemannSolver(const smf::PState &zL, const smf::PState &zR,
                  const StiffenedGas &eosL, const StiffenedGas &eosR,
                  double x_jump = 0.0);
    
    /// @brief Основной конструктор
    /// @param rL, uL, pL Плотность, скорость, давление слева
    /// @param rR, uR, pR Плотность, скорость, давление справа
    /// @param gL, p0L, e0L Параметра материала слева (двучленный УрС)
    /// @param gR, p0R, e0R Параметры материала справа (двучленный УрС)
    /// @param x_jump Положение разрыва
    RiemannSolver(
            double rL, double uL, double pL, double gL, double p0L, double e0L,
            double rR, double uR, double pR, double gR, double p0R, double e0R,
            double x_jump = 0.0);

    /// @brief Скорость звука от координаты и времени
    double sound_speed(double x, double t) const;

    /// @brief Плотность от координаты и времени
    double density(double x, double t) const;

    /// @brief Скорость от координаты и времени
    double velocity(double x, double t) const;

    /// @brief Давление от координаты и времени
    double pressure(double x, double t) const;

    /// @brief Энергия от координаты и времени
    double energy(double x, double t) const;

    /// @brief Характеристическая функция левой части области
    double fraction(double x, double t) const;

private:
    /// @brief Решение задачи Римана,
    /// инициализирует все поля класса
    void compute();

    // Следующие данные задаются в конструкторе

    // Положение начального разрыва
    double x_jump;

    // Параметры двучленного УрС
    double gL, p0L, e0L;
    double gR, p0R, e0R;

    // Начальные данные
    double rL, uL, pL;
    double rR, uR, pR;

    // Следующие поля инициализируются в compute()

    // Скорость звука
    double cL, cR;

    // Скорость и давление на контакте
    double U, P;

    // Параметры решения, плотность вокруг контактного разрыва
    double rl, rr;

    // Конфигурация решения (скорости УВ или ВР)
    double DL1, DL2, DR1, DR2;
};

} // namespace zephyr::math