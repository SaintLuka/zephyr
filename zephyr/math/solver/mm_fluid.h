#pragma once

#include <cmath>
#include <mutex>
#include <zephyr/mesh/mesh.h>
#include <zephyr/math/cfd/limiter.h>
#include <zephyr/phys/eos/eos.h>
#include <boost/format.hpp>
#include <zephyr/phys/tests/test_1D.h>
#include <zephyr/math/solver/riemann.h>
#include <zephyr/math/cfd/fluxes.h>

namespace zephyr::math {

using zephyr::mesh::Cell;
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
        std::array<double, mmf::Fractions::max_size> densities, speeds;

        [[nodiscard]] bool is_bad1() const {
            return get_pstate().is_bad();
        }

        [[nodiscard]] bool is_bad2() const {
            return half.is_bad();
        }

        [[nodiscard]] bool is_bad3() const {
            return next.is_bad();
        }

        [[nodiscard]] bool is_bad() const {
            return is_bad1() || is_bad3();
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
            densities.fill(0);
            speeds.fill(0);
        }
    };

    friend std::ostream &operator<<(std::ostream &os, const State &state) {
        os << boost::format(
                "Main state: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, temperature: %6%, energy: %7%, mass_frac: %8%\n") %
              state.rho % state.v.x() % state.v.y() % state.v.z() % state.p %
              state.t % state.e % state.mass_frac;
        os << boost::format(
                "Half state: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, temperature: %6%, energy: %7%, mass_frac: %8%\n") %
              state.half.density % state.half.velocity.x() % state.half.velocity.y() % state.half.velocity.z() %
              state.half.pressure % state.half.temperature % state.half.energy % state.half.mass_frac;
        os << boost::format(
                "Next state: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, temperature: %6%, energy: %7%, mass_frac: %8%\n") %
              state.next.density % state.next.velocity.x() % state.next.velocity.y() % state.next.velocity.z() %
              state.next.pressure % state.next.temperature % state.next.energy % state.next.mass_frac;
        return os;
    }

    /// @brief Получить экземпляр расширенного вектора состояния
    [[nodiscard]] static State datatype();

    /// @brief Конструктор класса
    explicit MmFluid(const phys::Materials &mixture, Fluxes flux = Fluxes::GODUNOV, double g = 0.0);

    /// @brief Число Куранта
    [[nodiscard]] double CFL() const;

    /// @brief Установить число Куранта
    void set_CFL(double CFL);

    void set_acc(int acc);

    void set_dim(int dim_);

    [[nodiscard]] double get_time() const;

    [[nodiscard]] size_t get_step() const;

    [[nodiscard]] std::string get_flux_name() const;

    /// @brief Шаг интегрирования на предыдущем вызове update()
    [[nodiscard]] double dt() const;

    void update(Mesh &mesh);

    /// @brief Установить флаги адаптации
    void set_flags(Mesh &mesh);

    /// @brief Распределитель данных при адаптации
    Distributor distributor() const;

private:

    void check_state(const mmf::PState &next, const mmf::QState &qc, const mmf::PState &old, const std::string &func_name, size_t cell_idx);

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    double compute_dt(Mesh &mesh);

    void compute_components_chars(Mesh &mesh);

    /// @brief Расчёт потоков
    void fluxes(Mesh &mesh);

    /// @brief Обновление ячеек
    void swap(Mesh &mesh);

    void compute_grad(Mesh &mesh, const std::function<mmf::PState(Cell &)> &to_state);

    void fluxes_stage1(Mesh &mesh);

    void fluxes_stage2(Mesh &mesh);

    mmf::Flux calc_flux_extra(Cell &state, bool from_begin);

public:

    ~MmFluid() = default;

protected:
    const phys::Materials mixture;
    int dim = 3;
    NumFlux::Ptr m_nf; ///< Метод расчёта потока
    double g = 0.0; ///< ускорение свободного падения, направлено против oy
    int m_acc = 1;
    double m_time = 0.0; ///< Прошедшее время
    size_t m_step = 0; ///< Количество шагов расчёта
    double m_CFL; ///< Число Куранта
    double m_dt; ///< Шаг интегрирования

    std::mutex out_mu{};
};

} // namespace zephyr