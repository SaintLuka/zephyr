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
        double rho; ///< плотность
        Vector3d v; ///< скорость
        double p; ///< давление
        double e; ///< энергия
        double t; ///< температура
        mmf::Fractions mass_frac; ///< доли веществ
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
            return mmf::PState(rho, v, p, e, t, mass_frac);
        }

        void set_state(const mmf::PState &pstate) {
            rho = pstate.density;
            v = pstate.velocity;
            p = pstate.pressure;
            e = pstate.energy;
            t = pstate.temperature;
            mass_frac = pstate.mass_frac;
        }
    };

    friend std::ostream &operator<<(std::ostream &os, const State &state) {
        os << boost::format(
                "State1: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, temperature: %6%, energy: %7%, mass_frac: %8%\n") %
              state.rho % state.v.x() % state.v.y() % state.v.z() % state.p %
              state.t % state.e % state.mass_frac;
        os << boost::format(
                "State2: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, temperature: %6%, energy: %7%, mass_frac: %8%\n") %
              state.next.density % state.next.velocity.x() % state.next.velocity.y() % state.next.velocity.z() %
              state.next.pressure % state.next.temperature % state.next.energy % state.next.mass_frac;
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
    void fluxes(Mesh &mesh) const;

    /// @brief Обновление ячеек
    void swap(Mesh &mesh);

    void compute_grad(Mesh &mesh,  const std::function<mmf::PState(zephyr::mesh::ICell &)> &to_state) const;

    void fluxes_stage1(Mesh &mesh) const;

    void fluxes_stage2(Mesh &mesh) const;

    mmf::Flux calc_flux_extra(ICell &cell, bool from_begin) const;

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