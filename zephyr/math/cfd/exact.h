#pragma once

#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/eos/eos.h>

namespace zephyr { namespace math {

/// @brief Класс с набором функций для получения точного решения 
/// задачи Римана о распаде разрыва. 
/// Наиболее эффективные (быстрые) функции являются статичными,
/// экземпляр класса требуется создавать в том случае, если необходимо
/// получить полное решение с зависимостью от времени.
class RiemannSolver {
public:
    using ref  = double&;
    using cref = const double&;

    /// @brief Начальная оценка давления на контакте
    /// @param rL, uL, pL Плотность, скорость, давление слева
    /// @param rR, uR, pR Плотность, скорость, давление справа
    /// @param cL, cR Скорость звука слева и справа
    /// @param p0L, p0R Параметры УрС
    /// @return Оценка давления на контакте
    static double contact_p_init(
            cref rL, cref uL, cref pL, cref cL, cref p0L,
            cref rR, cref uR, cref pR, cref cR, cref p0R);

    /// @brief Вычислить давление на контакте, используется
    /// итерационная процедура Годунова или аналог.
    /// @param rL, uL, pL Плотность, скорость, давление слева
    /// @param rR, uR, pR Плотность, скорость, давление справа
    /// @param gL, p0L Параметра материала слева (двучленный УрС)
    /// @param gR, p0R Параметры материала справа (двучленный УрС)
    /// @return Давление на контакте
    static double contact_p(
            cref rL, cref uL, cref pL, cref gL, cref p0L,
            cref rR, cref uR, cref pR, cref gR, cref p0R);

    /// @brief Вычислить давление на контакте, используется
    /// итерационная процедура Годунова или аналог. Отличается от
    /// предыдущей функции тем, что использует известные скорости
    /// звука слева и справа.
    /// @param cL, cR Скорость звука слева и справа, используются
    /// для начального приближения давления
    /// @param rL, uL, pL Плотность, скорость, давление слева
    /// @param rR, uR, pR Плотность, скорость, давление справа
    /// @param gL, p0L Параметра материала слева (двучленный УрС)
    /// @param gR, p0R Параметры материала справа (двучленный УрС)
    /// @return Давление на контакте
    static double contact_p(
            cref rL, cref uL, cref pL, cref cL, cref gL, cref p0L,
            cref rR, cref uR, cref pR, cref cR, cref gR, cref p0R);

    /// @brief Вычислить скорость на контакте
    /// @param uL, pL, aL Скорость, давление, массовая скорость слева
    /// @param rR, uR, pR Скорость, давление, массовая скорость справа
    /// @return Скорость на контакте
    static double contact_u(
            cref uL, cref pL, cref aL,
            cref uR, cref pR, cref aR);

    /// @brief Одноматериальный конструктор
    /// @param zL Вектор состояния слева
    /// @param zR Вектор состояния справа
    /// @param eos Уравнение состояния
    RiemannSolver(const smf::PState &zL, const smf::PState &zR,
                  const phys::Eos &eos, double x_jump = 0.0);

    /// @brief Двухматериальный конструктор
    /// @param zL Вектор состояния слева
    /// @param zR Вектор состояния справа
    /// @param eosL Уравнение состояния слева
    /// @param eosR Уравнение состояния справа
    RiemannSolver(const smf::PState &zL, const smf::PState &zR,
                  const phys::Eos &eosL, const phys::Eos &eosR,
                  double x_jump = 0.0);
    
    /// @brief Основной конструктор
    /// @param rhoL, uL, pL Плотность, скорость, давление слева
    /// @param rhoR, uR, pR Плотность, скорость, давление справа
    /// @param gL, p0L, e0L Параметра материала слева (двучленный УрС)
    /// @param gR, p0R, e0R Параметры материала справа (двучленный УрС)
    /// @return Давление на контакте
    RiemannSolver(
            double rhoL, double uL, double pL, double gL, double p0L, double e0L,
            double rhoR, double uR, double pR, double gR, double p0R, double e0R,
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

private:
    /// @brief Решение задачи Римана
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

    // Параметры решения, значения вокруг контактного разрыва
    double rl, rr, cl, cr;

    // Конфигурация решения (характеристики)
    double DL1, DL2, DR1, DR2;
};

}
}