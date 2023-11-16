#pragma once

#include <cmath>
#include <zephyr/mesh/mesh.h>
#include <zephyr/math/cfd/limiter.h>
#include <zephyr/phys/eos/eos.h>
#include <boost/format.hpp>
#include <zephyr/phys/tests/classic_test.h>
#include <zephyr/math/solver/riemann.h>
#include <zephyr/math/cfd/fluxes.h>

namespace zephyr::math {

using zephyr::mesh::ICell;
using zephyr::mesh::Mesh;
using zephyr::mesh::Distributor;
using zephyr::geom::Vector3d;

class MmFluid {
public:

    /// @brief Расширенный вектор состояния на котором решается задача
    struct State {
        double rho1; ///< плотность
        Vector3d v1; ///< скорость
        double p1; ///< давление
        double e1; ///< энергия
        double t1; ///< температура
        mmf::Fractions mass_frac1; ///< доли веществ
        mmf::PState half, next;
        mmf::PState d_dx, d_dy, d_dz;

        [[nodiscard]] bool is_bad2() const {
            return next.is_bad();
        }

        [[nodiscard]] bool is_bad1() const {
            return get_pstate().is_bad();
        }

        [[nodiscard]] bool is_bad() const {
            return is_bad1() || is_bad2();
        }

        [[nodiscard]] mmf::PState get_pstate() const {
            return mmf::PState(rho1, v1, p1, e1, t1, mass_frac1);
        }

        void set_state(const mmf::PState &pstate) {
            rho1 = pstate.density;
            v1 = pstate.velocity;
            p1 = pstate.pressure;
            e1 = pstate.energy;
            t1 = pstate.temperature;
            mass_frac1 = pstate.mass_frac;
        }
    };

    friend std::ostream &operator<<(std::ostream &os, const State &state) {
        os << boost::format(
                "State: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, temperature: %6%, energy: %7%, mass_frac: %8%\n") %
              state.rho1 % state.v1.x() % state.v1.y() % state.v1.z() % state.p1 %
              state.t1 % state.e1 % state.mass_frac1;
        return os;
    }

    /// @brief Получить экземпляр расширенного вектора состояния
    [[nodiscard]] static State datatype();

    /// @brief Конструктор класса
    explicit MmFluid(const phys::Materials &mixture, Fluxes flux);

    /// @brief Число Куранта
    [[nodiscard]] double CFL() const;

    /// @brief Установить число Куранта
    void set_CFL(double CFL);

    void set_acc(int acc);

    [[nodiscard]] double get_time() const;

    [[nodiscard]] size_t get_step() const;

    [[nodiscard]] std::string get_flux_name() const;

    /// @brief Шаг интегрирования на предыдущем вызове update()
    [[nodiscard]] double dt() const;

    void update(Mesh &mesh);

private:
    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    double compute_dt(Mesh &mesh);

    /// @brief Расчёт потоков
    void fluxes(Mesh &mesh);

    /// @brief Обновление ячеек
    void swap(Mesh &mesh);

    void compute_grad(Mesh &mesh);

    void fluxes_stage1(Mesh &mesh);

    void fluxes_stage2(Mesh &mesh);

    mmf::Flux calc_flux_extra(ICell &cell);

public:

    ~MmFluid() = default;

protected:
    const phys::Materials mixture;
    NumFlux::Ptr m_nf; ///< Метод расчёта потока
    int m_acc = 1;
    double m_time = 0.0; ///< Прошедшее время
    size_t m_step = 0; ///< Количество шагов расчёта
    double m_CFL; ///< Число Куранта
    double m_dt; ///< Шаг интегрирования
};

} // namespace zephyr