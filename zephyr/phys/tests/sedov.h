#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/matter/eos/ideal_gas.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;

/// Автомодельное решение задачи о точечном взрыве в сферических координатах.
/// Аналитическое решение из следующих источников:
/// https://github.com/brantr/sedov_taylor
/// https://en.wikipedia.org/wiki/Taylor%E2%80%93von_Neumann%E2%80%93Sedov_blast_wave
/// Но многие обозначения изменены, узнается с трудом.
class Sedov3D {
public:
    struct params {
        double gamma = 1.4;
        double rho0  = 1.0;
        double E     = 1.0;
    };
    
    /// @brief Конструктор
    Sedov3D(const params& p);

    /// @brief Радиус ударной волны
    double r_shock(double t) const;

    /// @brief Время, к которому УВ достигает радиус r
    double time_by_radius(double r) const;

    /// @brief Квадрат скорости звука
    double c_squared(double r, double t) const;

    /// @brief Плотность
    double density(double r, double t) const;

    /// @brief Модуль скорости
    double velocity(double r, double t) const;

    /// @brief Давление
    double pressure(double r, double t) const;

    /// @brief Внутренняя энергия
    double energy(double r, double t) const;

protected:
    double gamma;  ///< Показатель адиабаты
    double rho0;   ///< Плотность в невозмущенной области
    double E;      ///< Энерговыделение

    double beta;   ///< Параметр решения

    /// @brief Узлы сплайнов
    std::vector<double> m_xi;
    std::vector<double> m_Vs;
    std::vector<double> m_Gs;
    std::vector<double> m_Zs;
    std::vector<double> m_Ps; // xi^2 * G(xi) * Z(xi)

    double V_xi(double xi) const;
    double Z_xi(double xi) const;
    double G_xi(double xi) const;
    double P_xi(double xi) const;
};

} // namespace zephyr::phys