#pragma once

#include <zephyr/mesh/mesh.h>
#include <zephyr/math/cfd/limiter.h>
#include "zephyr/phys/eos/eos.h"
#include "zephyr/phys/tests/classic_test.h"
#include "riemann.h"
#include <zephyr/math/cfd/fluxes.h>

namespace zephyr { namespace math {

using zephyr::mesh::EuCell;
using zephyr::mesh::EuMesh;
using zephyr::mesh::Distributor;
using zephyr::geom::Vector3d;

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;


class SmFluid {
public:

    struct State {
        double rho;
        Vector3d v;
        double p;
        double e;
        PState half, next;
        PState d_dx, d_dy, d_dz;

        PState get_state() const {
            return {rho, v, p, e};
        }

        void set_state(const PState& z) {
            rho = z.density;
            v = z.velocity;
            p = z.pressure;
            e = z.energy;
        }     
    };

    /// @brief Получить экземпляр расширенного вектора состояния
    [[nodiscard]] static State datatype();

    ///@brief Конструктор класса, параметры по умолчанию
    SmFluid(const phys::Eos &eos, Fluxes flux);

    ///@brief
    void init_cells(EuMesh &mesh, const phys::ClassicTest &test);

    ///@brief 
    double CFL() const;

    ///@brief
    std::string get_flux_name() const;

    ///@brief
    void set_CFL(double CFL);

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    void compute_dt(EuMesh &mesh);

    /// @brief Один шаг интегрирования по времени
    void update(EuMesh& mesh);

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

    // Дифференциальный поток
    Flux calc_flux(EuCell &cell);

    // С экстраполяцией
    Flux calc_flux_extra(EuCell &cell, bool from_begin); 

    ///@brief вычисление градиента
    void compute_grad(EuMesh &mesh,  const std::function<smf::PState(zephyr::mesh::EuCell &)> &to_state) const;

    /// @brief установить порядок метода
    void set_accuracy(int acc);

    ///@brief Стадия 1
    void fluxes_stage1(EuMesh &mesh);

    ///@brief Стадия 2
    void fluxes_stage2(EuMesh &mesh);

    /// @brief Распределитель данных при адаптации
    Distributor distributor() const;

protected:
    const phys::Eos &m_eos;
    NumFlux::Ptr m_nf; ///< Метод расчёта потока
    double m_time = 0.0; ///< Прошедшее время
    size_t m_step = 0; ///< Количество шагов расчёта
    double m_CFL; ///< Число Куранта
    double m_dt; ///< Шаг интегрирования
    int m_acc = 2; ///< Порядок
};
}
}
