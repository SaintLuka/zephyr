#pragma once

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/math/cfd/limiter.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/eos/ideal_gas.h>

namespace zephyr { namespace math {

using zephyr::mesh::EuCell;
using zephyr::mesh::EuMesh;
using zephyr::mesh::Distributor;
using zephyr::geom::Vector3d;

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

/// @brief GK_solver 2D???

class SmFluid {
public:

    struct State {
        double rho1, rho2;
        Vector3d v1, v2;
        double p1, p2;
        double e1, e2;

        PState get_state1() const {
            return {rho1, v1, p1, e1};
        }

        void set_state2(const PState& z) {
            rho2 = z.density;
            v2 = z.velocity;
            p2 = z.pressure;
            e2 = z.energy;
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
    SmFluid();

    void compute_grad(EuCell &cell, int stage);

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    double compute_dt(EuCell& cell);

    /// @brief Один шаг интегрирования по времени
    void update(EuMesh& mesh, IdealGas &eos);

    /// @brief Векторное поле скорости
    /// @details Виртуальная функция, следует унаследоваться от класса
    /// Convection и написать собственную функцию скорости
    virtual Vector3d velocity(const Vector3d& c) const;

protected:

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double dt() const;

    /// @brief 
    void fluxes(EuCell& cell, int stage);

    void solution_step();

public:
    double m_dt;        ///< Шаг интегрирования

protected:
    
    //double m_dt;        ///< Шаг интегрирования
    Limiter m_limiter;  ///< Ограничитель
    double m_CFL;
};
}
}

