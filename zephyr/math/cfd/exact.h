#pragma once

#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/eos/stiffened_gas.h>

namespace zephyr { namespace math {

/// @brief Класс с точным решением задачи Римана.
/// Материалы слева и справа в общем случае могут быть различными.
class ExactRiemannSolver {
public:
    /// @brief Конструктор по умолчанию
    ExactRiemannSolver() = default;

    /// @brief Одноматериальный конструктор
    ExactRiemannSolver(const PState &zL, const PState &zR,
                       const phys::StiffenedGas &gas);

    /// @brief Двухматериальный конструктор
    ExactRiemannSolver(const PState &zL, const PState &zR,
                       const phys::StiffenedGas &gasL,
                       const phys::StiffenedGas &gasR);

    /// @brief Установить состояние слева
    void set_left(double density, Vector3d velocity, double pressure,
                  double energy, const phys::StiffenedGas &gas);

    /// @brief Установить состояние справа
    void set_right(double density, Vector3d velocity, double pressure,
                   double energy, const phys::StiffenedGas &gas);


    /// @brief Решить задачу Римана о распаде разрыва
    void compute();


    /// @brief Материал от координаты и времени
    int material(double x, double t) const;

    /// @brief Плотность от координаты и времени
    double density(double x, double t) const;

    /// @brief Скорость от координаты и времени
    Vector3d velocity(double x, double t) const;

    /// @brief Давление от координаты и времени
    double pressure(double x, double t) const;

    /// @brief Энергия от координаты и времени
    double energy(double x, double t) const;

private:
    // Начальные данные
    phys::StiffenedGas gasL, gasR;
    double rL, uL, vL, wL, eL, pL;
    double rR, uR, vR, wR, eR, pR;

    // Параметры решения, значения вокруг контактного разрыва
    double rl, rr, el, er, u, p;

    // Конфигурация решения
    double DL, DR; ///<
    double cL, cR;
    double cl, cr;

};

}
}