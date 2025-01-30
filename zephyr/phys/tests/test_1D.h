/// @file Набор одномерных газодинамических тестов.
/// Что есть сейчас:
///   SodTest       -- Классический тест Сода
///   ToroTest      -- Тесты Торо из монографии
///   ShockWave     -- Простая ударная волна, задается числом Маха
///   ShuOsher      -- Тест Шу-Ошера, волна падает на возмущение
///   RarefiedWater -- Тест с разреженной водой и газом
///   Multimat1D    -- Многоматериальные одномерные тесты

#pragma once

#include <zephyr/phys/tests/ivp.h>
#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/math/solver/riemann.h>

namespace zephyr::phys {

using zephyr::math::RiemannSolver;

/// @class Абстрактный класс одномерного теста
class Test1D : public IVP {
public:
    /// @brief Левая граница области
    virtual double xmin() const = 0;

    /// @brief Правая граница области
    virtual double xmax() const = 0;

    /// @brief Получить название теста
    std::string name() const override {
        return "Test1D";
    }
};

/// @class Классический одномерный тест Сода
class SodTest : public Test1D {
public:
    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия

    /// @brief Точное решение задачи Римана
    RiemannSolver exact;

    /// @brief Конструктор
    SodTest();

    std::string name() const final { return "SodTest"; }

    /// @brief Левая граница области
    double xmin() const final { return 0.0; }

    /// @brief Правая граница области
    double xmax() const final { return 1.0; }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    // Начальные данные

    double density(const Vector3d &r) const final;

    Vector3d velocity(const Vector3d &r) const final;

    double pressure(const Vector3d &r) const final;

    // Точное решение

    double density_t(const Vector3d &r, double t) const final;

    Vector3d velocity_t(const Vector3d &r, double t) const final;

    double pressure_t(const Vector3d &r, double t) const final;
};

/// @class Набор тестов на распад разрыва из монографии Торо
/// (глава 10 и 4.3.3 Numerical Tests)
/// E.F. Toro. Riemann Solvers and Numerical Methods for Fluid Dynamics.
class ToroTest : public Test1D {
protected:
    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия

    /// @brief Точное решение задачи Римана
    RiemannSolver exact;

    // Обновляет энергии и точное решение после изменения
    // параметров или уравнений состояния
    void update();

public:
    /// @brief Конструктор
    /// @param num Номер теста 1..7
    explicit ToroTest(int num, bool multimat = false);

    /// @brief Выравнять температуру на разрыве путем
    /// настройки удельной теплоемкости
    void adjust_cv();

    std::string name() const final { return "ToroTest";}

    /// @brief Левая граница области
    double xmin() const final { return 0.0; }

    /// @brief Правая граница области
    double xmax() const final { return 1.0; }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    // Начальные данные

    int index(const Vector3d& r) const final;

    double density(const Vector3d& r) const final;

    Vector3d velocity(const Vector3d& r) const final;

    double pressure(const Vector3d& r) const final;

    // Точное решение

    int index_t(const Vector3d& r, double t) const final;

    double density_t(const Vector3d &r, double t) const final;

    Vector3d velocity_t(const Vector3d &r, double t) const final;

    double pressure_t(const Vector3d &r, double t) const final;
};

/// @brief Простая ударная волна, задается числом Маха и параметрами
/// перед фронтом ударной волны
class ShockWave : public Test1D {
public:
    double x_jump;  ///< Положение разрыва
    double length;  ///< Длина области
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление

    double speed;   ///< Скорость разрыва

    // Параметры перед фронтом УВ
    struct Params {
        double rho    = 1.0;
        double P      = 1.0;
        double gamma  = 1.4;
    };

    /// @brief Конструктор
    /// @param Ms Число Маха (при Ms < 0.0, направлена влево)
    /// @param x_jump Начальное положение разрыва
    /// @param L длина области
    explicit ShockWave(double Ms = 6, double x_jump = 0.1, double L = 1.0,
                       Params params = Params{.rho=1.0, .P=1.0, .gamma=1.4});

    std::string name() const final { return "Shock Wave";}

    /// @brief Левая граница области
    double xmin() const final { return 0.0; }

    /// @brief Правая граница области
    double xmax() const final { return length; }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }


    // Начальные данные

    double density(const Vector3d& r) const final;

    Vector3d velocity(const Vector3d& r) const final;

    double pressure(const Vector3d& r) const final;

    // Точное решение

    double density_t(const Vector3d &r, double t) const final;

    Vector3d velocity_t(const Vector3d &r, double t) const final;

    double pressure_t(const Vector3d &r, double t) const final;
};

/// @brief Тест Шу-Ошера
/// C.-W. Shu and S. Osher. Efficient Implementation of Essentially
/// Non-oscillatory Shock-Capturing Schemes, II (1988)
class ShuOsherTest : public Test1D {
public:
    double epsilon;
    double rL, uL, pL;
    double rR, uR, pR;

    const double x_jump = -4.0;
    double finish = 1.8;

    /// @brief Конструктор
    ShuOsherTest();

    /// @brief Получить название теста
    std::string name() const final { return "Shu-Osher Test"; };

    /// @brief Левая граница области
    double xmin() const final { return -5.0; }

    /// @brief Правая граница области
    double xmax() const final { return +5.0; }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    // Начальные данные

    double density(const Vector3d& r) const final;

    Vector3d velocity(const Vector3d& r) const final;

    double pressure(const Vector3d& r) const final;
};

/// @brief Тест с разреженной водой.
/// Отрицательное начальное давление воды, контакт с газом.
class RarefiedWater : public Test1D {
public:
    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия

    /// @brief Конструктор
    RarefiedWater();

    std::string name() const final { return "Rarefied water";}

    /// @brief Левая граница области
    double xmin() const final { return 0.0; }

    /// @brief Правая граница области
    double xmax() const final { return 1.0_cm; }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    // Начальные данные

    int index(const Vector3d& r) const final;

    double density(const Vector3d &r) const final;

    Vector3d velocity(const Vector3d &r) const final;

    double pressure(const Vector3d &r) const final;
};

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

    // Начальные данные

    double density(const Vector3d &r) const final;

    Vector3d velocity(const Vector3d &r) const final;

    double pressure(const Vector3d &r) const final;

    double energy(const Vector3d &r) const final;

    double temperature(const Vector3d &r) const final;
};

} // namespace zephyr::phys