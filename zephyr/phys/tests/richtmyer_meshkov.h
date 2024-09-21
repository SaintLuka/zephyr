#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/geom/primitives/boundary.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;
using zephyr::geom::Boundary;
using zephyr::geom::generator::Rectangle;

class RichtmyerMeshkov {
public:
    Rectangle generator; ///< Сеточный генератор

    IdealGas::Ptr eosL;  ///< Используемый УрС
    IdealGas::Ptr eosR;  ///< Используемый УрС

    double finish;       ///< Конечный момент времени
    double x_contact;    ///< Положение контакта
    double x_jump;       ///< Положение ударной волны
    double y_width;      ///< Толщина трубки
    double amplitude;    ///< Амплитуда возмущения
    int    n_peaks;      ///< Число возмущений

    double pL, rL, eL, uL;  ///< За фронтом УВ
    double pR, rR, eR, uR;  ///< Перед фронтом УВ
    double p0, r0, e0, u0;  ///< Область за контактом


    /// @brief Конструктор
    /// @param Ms Число Маха
    explicit RichtmyerMeshkov(double Ms = 3.0) {
        double gamma = 1.4;
        eosL = IdealGas::create(gamma);
        eosR = eosL;

        // Невозмущенная область
        p0 = 1.0;
        r0 = 0.01;
        u0 = 0.0;
        e0 = eosR->energy_rp(r0, p0);

        // Перед фронтом УВ
        pR = 1.0;
        rR = 1.0;
        uR = 0.0;
        eR = eosR->energy_rp(rR, pR);

        // За фронтом УВ
        pL = pR * (2 * gamma * Ms * Ms - gamma + 1) / (gamma + 1);
        rL = rR * (gamma + 1) * Ms * Ms / (2 + (gamma - 1) * Ms * Ms);
        uL = 2.0 / Ms * std::sqrt(gamma * pR / rR) * (Ms * Ms - 1) / (gamma + 1);
        eL = eosL->energy_rp(rL, pL);

        x_jump = 0.5;
        x_contact = 1.0;
        amplitude = 0.02;
        y_width = 1.0;
        n_peaks = 3;

        finish = 3.0;

        generator = Rectangle(0.0, 10.0, -0.5 * y_width, 0.5 * y_width);
        generator.set_boundaries({.left=Boundary::ZOE, .right=Boundary::WALL,
                                  .bottom=Boundary::WALL, .top=Boundary::WALL});
    };

    /// @brief Название теста
    std::string name() const { return "Richtmyer-Meshkov instability";}

    /// @brief Конечный момент времени
    double max_time() const { return finish; }

    /// @brief Номер области
    /// 0 -- за фронтом ударной волны
    /// 1 -- между УВ и контактом
    /// 2 -- Невозмущенная область справа от контакта
    int region(const Vector3d& r) const {
        if (r.x() < x_jump) {
            return 0;
        }
        // длина волны возмущения
        double L = y_width / n_peaks;
        // s in [0, 1]
        double s = r.y() / L - std::floor(r.y() / L);
        // пила
        double x = x_contact + amplitude * (std::abs(s - 0.5) - 0.5);
        return r.x() < x ? 1 : 2;
    }

    /// @brief Начальная плотность
    double density(const Vector3d &r) const {
        switch (region(r)) {
            case 0:  return rL;
            case 1:  return rR;
            default: return r0;
        }
    };

    /// @brief Начальная скорость
    Vector3d velocity(const Vector3d &r) const {
        switch (region(r)) {
            case 0:  return {uL, 0.0, 0.0};
            case 1:  return {uR, 0.0, 0.0};
            default: return {u0, 0.0, 0.0};
        }
    };

    /// @brief Начальное давление
    double pressure(const Vector3d &r) const {
        switch (region(r)) {
            case 0:  return pL;
            case 1:  return pR;
            default: return p0;
        }
    };

    /// @brief Начальная внутренняя энергия
    double energy(const Vector3d &r) const {
        switch (region(r)) {
            case 0:  return eL;
            case 1:  return eR;
            default: return e0;
        }
    };

    ///@brief Получить используемый УрС
    IdealGas::Ptr get_eos(const Vector3d& r) const {
        switch (region(r)) {
            case 0:  return eosL;
            case 1:  return eosL;
            default: return eosR;
        }
    }
};

} // namespace zephyr::phys