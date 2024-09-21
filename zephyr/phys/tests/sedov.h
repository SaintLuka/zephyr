#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/phys/tests/test_1D.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;

// The Sedov (or Sedov-Taylor) blast wave is a standard hydrodynamics test problem
// A large amount of energy is placed into a very small volume, driving a spherical
// (or cylindrical in 2-d Cartesian coordinates) blast wave
class SedovBlast {
public:
    int constructor_type;

    IdealGas::Ptr eos;       ///< Используемый УрС
    double r_jump;           ///< Положение разрыва
    double finish;           ///< Конечный момент времени
    double r_init, r_amb;    ///< Плотность
    Vector3d v_init, v_amb;  ///< Скорость
    double p_init, p_amb;    ///< Давление
    double E;                ///< Начальная энергия
    double V_init;           ///< Начальный объем
    Rectangle generator;     ///< Сеточный генератор


    /// @brief Конструктор
    SedovBlast(double gamma=1.4) : constructor_type(1) {
        eos = IdealGas::create(gamma);

        E = 1; // у.е.
        r_jump = 0.02;
        V_init = M_PI * r_jump * r_jump;
        
        p_amb = 0.06;
        p_init = eos->pressure_re (1.0 / V_init, E);
        
        r_amb = 1.0;
        r_init = 1.0;

        v_amb  = Vector3d::Zero();
        v_init = Vector3d::Zero();

        finish = 0.15;

        generator = Rectangle(-0.5, 0.5, -0.5, 0.5);
        generator.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
                              .bottom=Boundary::ZOE, .top=Boundary::ZOE});
    }

    /// @brief Задать число Маха и начальный радиус ударной волны
    SedovBlast(double Ms, double x_jump, double gamma=1.4) : r_jump(x_jump), constructor_type(2) {
        eos = IdealGas::create(gamma);

        r_amb = 1.225;
        p_amb = 1;
        v_amb = {0,0,0};

        r_init = r_amb * (eos->gamma + 1) * Ms * Ms / (2 + (eos->gamma - 1) * Ms * Ms);
        p_init = p_amb * (1 + 2 * eos->gamma * (Ms * Ms - 1) / (eos->gamma + 1));
        double radV = Ms * sqrt(eos->gamma * p_amb / r_amb) * (1 - r_amb / r_init);
        v_init = {radV / sqrt(2), radV / sqrt(2), 0};

        finish = 0.5;

        generator = Rectangle(-0.25, 0.25, -0.25, 0.25);
        generator.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
                                  .bottom=Boundary::ZOE, .top=Boundary::ZOE});
    }

    std::string name() const { return "SedovBlast";}

    /// @brief Конечный момент времени
    double max_time() const { return finish; }

    ///@brief Получить положение разрыва
    double get_x_jump() const { return r_jump; }

    /// @brief Начальная плотность
    double density(const Vector3d &r) const {
        if (constructor_type < 2) {
            return r.norm() < r_jump ? r_init : r_amb;
        } else {
            return r.norm() < r_jump ? std::pow((r.norm() / r_jump), 3 / (eos->gamma - 1)) : r_amb;
        }
    }

    /// @brief Начальное давление
    double pressure(const Vector3d &r) const {
        if (constructor_type < 2) {
            return r.norm() < r_jump ? p_init : p_amb;
        } else {
            return r.norm() < r_jump ? p_init * std::pow(3, std::pow(r.norm() / r_jump, 5)) : p_amb;
        }
    }

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &r) const {
        if (constructor_type < 2) {
            return r.norm() < r_jump ? v_init : v_amb;
        } else {
            return r.norm() < r_jump ? v_init * (r.norm() / r_jump) : v_amb;
        }
    }


    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const { 
        return eos->energy_rp(density(r), pressure(r));
    }

    /// @brief Уравнение состояния
    Eos::Ptr get_eos(const Vector3d &r) const {
        return eos;
    }

    ~SedovBlast() = default;
};

} // namespace zephyr::phys