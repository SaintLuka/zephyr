#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;

/// "Автомодельное" решение в сферических координатах.
/// Аналитическое решение из следующих источников:
/// https://github.com/brantr/sedov_taylor
/// https://en.wikipedia.org/wiki/Taylor%E2%80%93von_Neumann%E2%80%93Sedov_blast_wave
class SedovBlast3D {
public:
    IdealGas::Ptr eos;  ///< Используемый УрС
    double finish;      ///< Конечный момент времени (радиус = 1)

    struct params {
        double gamma = 1.4;
        double rho0  = 1.0;
        double E     = 1.0;
    };

    /// @brief Конструктор
    SedovBlast3D(const params& p);

    std::string name() const { return "SedovBlast3D"; }

    /// @brief Конечный момент времени
    double max_time() const { return finish; }

    /// @brief Радиус ударной волны
    double r_shock(double t) const;

    /// @brief Время, к которому УВ достигает радиус r
    double time_by_radius(double r) const;

    /// @brief Плотность
    double density(const Vector3d &r, double t) const;

    /// @brief Давление
    double pressure(const Vector3d &r, double t) const;

    /// @brief Скорость
    Vector3d velocity(const Vector3d &r, double t) const;

    /// @brief Внутренняя энергия
    double energy(const Vector3d &r, double t) const;

protected:
    double gamma;  ///< Показатель адиабаты
    double rho0;   ///< Плотность в невозмущенной области
    double E;      ///< Энерговыделение

    /// @brief Параметры решения
    double v1, v2, v3, v4, v5;
    double beta;

    double eta_V(double V) const;
    double V_eta(double eta) const;

    double V_xi(double xi) const;
    double Z_xi(double xi, double V = NAN) const;
    double G_xi(double eta, double V = NAN) const;

    /// @brief Squared sound speed
    double c_squared(double r, double t) const;

    double rho(double r, double t) const;

    double v(double r, double t) const;

    double p(double r, double t) const;
};

} // namespace zephyr::phys