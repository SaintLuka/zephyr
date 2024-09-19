#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/tests/test_1D.h>
#include <zephyr/phys/tests/test_2D.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;
using zephyr::geom::generator::Rectangle;

/// @class Набор тестов на распад разрыва из монографии Торо (глава 10) и (4.3.3 Numerical Tests)
/// E.F. Toro. Riemann Solvers and Numerical Methods for Fluid Dynamics.
class ToroTest : public Test1D {
public:

    IdealGas::Ptr eos_L;  ///< Используемый УрС
    IdealGas::Ptr eos_R; ///< Используемый УрС

    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия


    /// @brief Конструктор
    /// @param num Номер теста 1..7
    explicit ToroTest(int num, bool multimat = false);

    /// @brief Симметрично отразить начальные условия
    void inverse();

    std::string get_name() const final { return "ToroTest";}

    /// @brief Левая граница области
    double xmin() const final { return 0.0; }

    /// @brief Правая граница области
    double xmax() const final { return 1.0; }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    ///@brief Получить используемый УрС
    Eos::Ptr get_eos() const final {
        return eos_L;
    }

    ///@brief Получить положение разрыва
    double get_x_jump() const final { return x_jump; }


    /// @brief Начальная плотность
    double density(double x) const final {
        return x < x_jump ? rL : rR;
    }

    /// @brief Начальная скорость
    Vector3d velocity(double x) const final {
        return { x < x_jump ? uL : uR, 0.0, 0.0};
    }

    /// @brief Начальное давление
    double pressure(double x) const final {
        return x < x_jump ? pL : pR;
    }

    /// @brief Начальная внутренняя энергия
    double energy(double x) const final {
        return x < x_jump ? eL : eR;
    }

    /// @brief Версия для многоматериальных
    Eos::Ptr get_eos(double x) const final {
        return x < x_jump ? eos_L : eos_R;
    }


    /// @brief Начальная плотность
    double density(const Vector3d& r) const final {
        return density(r.x());
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d& r) const final {
        return velocity(r.x());
    }

    /// @brief Начальное давление
    double pressure(const Vector3d& r) const final {
        return pressure(r.x());
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d& r) const final {
        return energy(r.x());
    }

    /// @brief Версия для многоматериальных
    Eos::Ptr get_eos(const Vector3d &r) const final {
        return get_eos(r.x());
    }

    ~ToroTest() = default;
};


/// @class Тесты Торо, повернутые на угол
class ToroTest2D {
public:
    ToroTest test1D;
    double angle;
    Rectangle m_gen;

    /// @brief Конструктор
    /// @param num Номер теста 1..7
    explicit ToroTest2D(int num, double angle, bool multimat = false)
        : test1D(num, multimat), angle(angle) {
        m_gen = Rectangle(-5.0, 5.0, -5.0, 5.0);
        m_gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
                              .bottom=Boundary::ZOE, .top=Boundary::ZOE});
    }

    Rectangle generator() const {
        return m_gen;
    }

    /// @brief Конечный момент времени
    double max_time() const { return 5.0 * test1D.finish; }

    ///@brief Получить используемый УрС
    const Eos& get_eos() const { return *test1D.eos_L; }

    ///@brief Получить положение разрыва
    double get_x_jump() const { return test1D.x_jump; }

    Vector3d rotate(const Vector3d& r) const {
        Vector3d res = {
                +std::cos(angle) * r.x() + std::sin(angle) * r.y(),
                -std::sin(angle) * r.x() + std::cos(angle) * r.y(),
                0.0
        };
        return res;
    }

    /// @brief Начальная плотность
    double density(const Vector3d &r) const {
        return test1D.density(rotate(r));
    }

    /// @brief Начальное давление
    double pressure(const Vector3d &r) const {
        return test1D.pressure(rotate(r));
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &r) const {
        return test1D.velocity(rotate(r));
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const {
        return test1D.energy(rotate(r));
    }

    ~ToroTest2D() = default;
};

} // namespace zephyr::phys
