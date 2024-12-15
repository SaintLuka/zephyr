#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/stiffened_gas.h>
#include <zephyr/phys/tests/test_1D.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;

/// @class Похоже на одномерные тесты Торо, но с многоматериальной постановкой.
/// x_jump -- положение разрыва, x_cont -- граница материалов.
/// Большое число одномерных тестов, часто с точным решением.
class Multimat1D : public Test1D {
public:
    StiffenedGas::Ptr eos_L;    ///< Используемый УрС
    StiffenedGas::Ptr eos_R;    ///< Используемый УрС

    double x_jump;              ///< Положение разрыва
    double x_cont;              ///< Граница материалов
    double finish;              ///< Конечный момент времени
    double rL1, rL2, rR1, rR2;  ///< Плотность
    double uL1, uL2, uR1, uR2;  ///< Скорость
    double pL1, pL2, pR1, pR2;  ///< Давление
    double eL1, eL2, eR1, eR2;  ///< Внутренняя энергия
    double tL1, tL2, tR1, tR2;  ///< Температура

    double width; ///< Ширина размазывания


    /// @brief Конструктор
    /// @param num Номер теста
    /// @param mat Конфигурация материалов
    /// @param c Вариант постановки (case)
    explicit Multimat1D(int num, int mat, int c);

    /// @brief Название теста
    std::string name() const final { return "Multimat1D"; }

    /// @brief Левая граница области
    double xmin() const final { return 0.0; }

    /// @brief Правая граница области
    double xmax() const final { return 1.0; }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    ///@brief Получить положение разрыва
    double get_x_jump() const final { return x_jump; }


    /// @brief Начальная плотность
    double density(const Vector3d &r) const final;

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &r) const final;

    /// @brief Начальная температура
    double pressure(const Vector3d &r) const final;

    /// @brief Начальное давление
    double temperature(const Vector3d &r) const;

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const final;

    /// @brief Уравнение состояния
    Eos::Ptr get_eos(const Vector3d &r) const final;

    /// @brief Доля материала
    double fraction(const Vector3d& r, int mat) const;

};

} // namespace zephyr::phys