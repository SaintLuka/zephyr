#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/tests/test_1D.h>

#include <zephyr/geom/generator/rectangle.h>

namespace zephyr::phys {

using namespace zephyr::geom;
using zephyr::geom::generator::Rectangle;

/// @brief Простая ударная волна
class ShockWave : public Test1D {
public:
    IdealGas::Ptr eos;   ///< Используемый УрС

    double x_jump;  ///< Положение разрыва
    double finish;  ///< Конечный момент времени
    double rL, rR;  ///< Плотность
    double uL, uR;  ///< Скорость
    double pL, pR;  ///< Давление
    double eL, eR;  ///< Внутренняя энергия

    /// @brief Конструктор
    ShockWave(double Ms=6, double x_jump=0.1, double finish=0.1)
        : x_jump(x_jump), finish(finish) {

        double gamma = 1.4;
        eos = IdealGas::create(gamma);

        pR = 1.0;
        rR = 1.0;
        uR = 0.0;
        
        pL = pR * (2 * gamma * Ms * Ms  - gamma + 1) / (gamma + 1) ;
        rL = rR * ( gamma + 1 ) * Ms * Ms / ( 2 + ( gamma - 1 ) * Ms * Ms );
        uL = 2 / Ms * std::sqrt( gamma * pR / rR ) * ( Ms * Ms - 1 ) / (gamma + 1 );

        eL = eos->energy_rp(rL, pL);
        eR = eos->energy_rp(rR, pR);
    }

    /// @brief Симметрично отразить начальные условия
    void inverse() {
        std::swap(rL, rR);
        std::swap(uL, uR);
        std::swap(pL, pR);
        std::swap(eL, eR);
        uL *= -1.0;
        uR *= -1.0;
    }

    std::string name() const final { return "Shock Wave";}

    /// @brief Левая граница области
    double xmin() const final { return 0.0; }

    /// @brief Правая граница области
    double xmax() const final { return 1.5; }

    /// @brief Конечный момент времени
    double max_time() const final { return finish; }

    ///@brief Получить положение разрыва
    double get_x_jump() const final { return x_jump; }


    /// @brief Начальная плотность
    double density(const Vector3d& r) const final {
        return r.x() < x_jump ? rL : rR;
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d& r) const final {
        return {r.x() < x_jump ? uL : uR, 0, 0};
    }

    /// @brief Начальное давление
    double pressure(const Vector3d& r) const final {
        return r.x() < x_jump ? pL : pR;
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d& r) const final {
        return r.x() < x_jump ? eL : eR;
    }

    /// @brief Уравнение состояния
    Eos::Ptr get_eos(const Vector3d& r) const final {
        return eos;
    }
};

/// @brief Косая ударная волна
class SkewShockWave {
public:
    ShockWave sw;
    double angle;         ///< Угол поворота
    Rectangle generator;  ///< Сеточный генератор

    /// @brief Конструктор
    SkewShockWave(double Ms=6, double angle = M_PI/6, double finish=0.2)
            : sw(Ms, 0.0, finish), angle(angle) {

        generator = Rectangle(-0.2, 2.3, 0.0, 1.0);
        generator.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
                                  .bottom=Boundary::WALL, .top=Boundary::ZOE});

    }

    std::string name() const { return "Skew Shock Wave";}

    /// @brief Левая граница области
    double xmin() const { return -0.2; }

    /// @brief Правая граница области
    double xmax() const { return 2.3; }

    /// @brief Конечный момент времени
    double max_time() const { return sw.finish; }

    /// @brief Вращение координат
    Vector3d rotate(const Vector3d& r, double phi) const {
        Vector3d res = {
                +std::cos(phi) * r.x() - std::sin(phi) * r.y(),
                +std::sin(phi) * r.x() + std::cos(phi) * r.y(),
                0.0
        };
        return res;
    }


    /// @brief Начальная плотность
    double density(const Vector3d& r) const {
        return sw.density(rotate(r, angle));
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d& r) const {
        return rotate(sw.velocity(rotate(r, angle)), -angle);
    }

    /// @brief Начальное давление
    double pressure(const Vector3d& r) const {
        return sw.pressure(rotate(r, angle));
    }

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d& r) const {
        return sw.energy(rotate(r, angle));
    }

    /// @brief Уравнение состояния
    Eos::Ptr get_eos(const Vector3d& r) const {
        return sw.get_eos(r);
    }
};

} // namespace zephyr::phys

