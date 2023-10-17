#pragma once

#include <zephyr/mesh/mesh.h>
#include <zephyr/math/cfd/limiter.h>
#include "zephyr/phys/eos/eos.h"
#include "zephyr/phys/tests/classic_test.h"
#include "riemann.h"
#include <zephyr/math/cfd/fluxes.h>

namespace zephyr { namespace math {

using zephyr::mesh::ICell;
using zephyr::mesh::Mesh;
using zephyr::mesh::Distributor;
using zephyr::geom::Vector3d;

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

/// @brief GK_solver 2D???

class SmFluid {
public:

    struct State {
        double rho1, rho2, rhoh, rhoh2;
        Vector3d v1, v2, vh, vh2;
        double p1, p2, ph, ph2;
        double e1, e2, eh, eh2;
        Flux P1, P2, P3, P4;

        // Производные
        double rhox; Vector3d vx; double px; double ex;
        double rhoy; Vector3d vy; double py; double ey;
        double rhoz; Vector3d vz; double pz; double ez; 

        PState zx() {
            return {rhox, vx, px, ex};
        }
        PState zy() {
            return {rhoy, vy, py, ey};
        }
        PState zz() {
            return {rhoz, vz, pz, ez};
        }

        PState set_zx(PState &z) {
            rhox = z.density;
            vx = z.velocity;
            px = z.pressure;
            ex = z.energy;
        }
        PState set_zy(PState &z) {
            rhoy = z.density;
            vy = z.velocity;
            py = z.pressure;
            ey = z.energy;
        }
        PState set_zz(PState &z) {
            rhoz = z.density;
            vz = z.velocity;
            pz = z.pressure;
            ez = z.energy;
        }

        PState get_state(int stage) const {
            if (stage == 0) {
                return {rho1, v1, p1, e1};
            }
            if (stage == 1) {
                return {rhoh, vh, ph, eh};
            }
            if (stage == 2) {
                return {rho2, v2, p2, e2};
            }
            if (stage == 3) {
                return {rhoh2, vh2, ph2, eh2};
            }
        }

        void set_state(const PState& z, int stage) {
            if (stage == 0) {
                rho1 = z.density;
                v1 = z.velocity;
                p1 = z.pressure;
                e1 = z.energy;
            }
            if (stage == 1) {
                rhoh = z.density;
                vh = z.velocity;
                ph = z.pressure;
                eh = z.energy;
            }
            if (stage == 2) {
                rho2 = z.density;
                v2 = z.velocity;
                p2 = z.pressure;
                e2 = z.energy;
            }
            if (stage == 3) {
                rhoh2 = z.density;
                vh2 = z.velocity;
                ph2 = z.pressure;
                eh2 = z.energy;
            }
        }

        void swap() {
            std::swap(rho1, rho2);
            std::swap(v1, v2);
            std::swap(p1, p2);
            std::swap(e1, e2);
        }
        
    };

    /// @brief Получить экземпляр расширенного вектора состояния
    static State datatype();

    ///@brief Конструктор класса, параметры по умолчанию
    SmFluid(const phys::Eos &eos, Fluxes flux);

    ///@brief
    void init_cells(Mesh &mesh, const phys::ClassicTest &test);

    ///@brief 
    double CFL() const;

    ///@brief
    std::string get_flux_name() const;

    ///@brief
    void set_CFL(double CFL);

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    double compute_dt(Mesh &mesh);

    /// @brief Один шаг интегрирования по времени
    void update(Mesh& mesh);

    /// @brief Векторное поле скорости
    /// @details Виртуальная функция, следует унаследоваться от класса
    /// Convection и написать собственную функцию скорости
    virtual Vector3d velocity(const Vector3d& c) const;

    /// @brief Шаг интегрирования на предыдущем вызове update()
    [[nodiscard]] double dt() const;

    ///@brief
    [[nodiscard]] double get_time() const;

    ///@brief
    [[nodiscard]] double get_step() const;

    Flux calc_flux(ICell &cell, int stage);

    void compute_grad(ICell &cell, int stage);

    ///@brief 
    void fluxes(Mesh &mesh);

    ///@brief
    void fluxes2(Mesh &mesh);

    void fluxes3(Mesh &mesh);

    ///@brief
    void fluxes4(Mesh &mesh);

protected:
    const phys::Eos &m_eos;
    NumFlux::Ptr m_nf; ///< Метод расчёта потока
    double m_time = 0.0; ///< Прошедшее время
    size_t m_step = 0; ///< Количество шагов расчёта
    double m_CFL; ///< Число Куранта
    double m_dt; ///< Шаг интегрирования
    Limiter m_limiter; ///< Лимитер
};
}
}

