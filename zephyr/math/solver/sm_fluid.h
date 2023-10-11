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
        double rho1, rho2, rhoh;
        Vector3d v1, v2, vh;
        double p1, p2, ph;
        double e1, e2, eh;

        PState get_state1() const {
            return {rho1, v1, p1, e1};
        };

        PState get_state_h() const {
            return {rhoh, vh, ph, eh};
        };

        void set_state2(const PState& z) {
            rho2 = z.density;
            v2 = z.velocity;
            p2 = z.pressure;
            e2 = z.energy;
        };

        void set_state_h(const PState &z) {
            rhoh = z.density;
            vh = z.velocity;
            ph = z.pressure;
            eh = z.energy;
        }

        void swap() {
            std::swap(rho1, rho2);
            std::swap(v1, v2);
            std::swap(p1, p2);
            std::swap(e1, e2);
        };
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

    ///@brief 
    void fluxes(Mesh &mesh);

    ///@brief
    void fluxes2(Mesh &mesh);

protected:
    const phys::Eos &m_eos;
    NumFlux::Ptr m_nf; ///< Метод расчёта потока
    double m_time = 0.0; ///< Прошедшее время
    size_t m_step = 0; ///< Количество шагов расчёта
    double m_CFL; ///< Число Куранта
    double m_dt; ///< Шаг интегрирования
};
}
}

