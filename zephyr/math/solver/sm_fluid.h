#pragma once

#include <cmath>
#include <zephyr/mesh/mesh.h>
#include <zephyr/math/cfd/limiter.h>
#include <zephyr/phys/eos/eos.h>
#include <boost/format.hpp>
#include <zephyr/math/cfd/fluxes.h>

namespace zephyr::math {

using zephyr::mesh::Cell;
using zephyr::mesh::Mesh;
using zephyr::mesh::Distributor;
using zephyr::geom::Vector3d;

class SmFluid {
public:

    /// @brief Расширенный вектор состояния на котором решается задача
    struct State {
        double rho; ///< плотность
        Vector3d v; ///< скорость
        double p; ///< давление
        double e; ///< энергия
        smf::PState half, next;
        smf::PState d_dx, d_dy, d_dz;

        [[nodiscard]] smf::PState get_pstate() const {
            return smf::PState(rho, v, p, e);
        }

        void set_state(const smf::PState &pstate) {
            rho = pstate.density;
            v = pstate.velocity;
            p = pstate.pressure;
            e = pstate.energy;
        }
    };

    friend std::ostream &operator<<(std::ostream &os, const State &state) {
        os << boost::format(
                "State1: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, energy: %6%\n") %
              state.rho % state.v.x() % state.v.y() % state.v.z() % state.p;
        os << boost::format(
                "State2: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, energy: %6%\n") %
              state.next.density % state.next.velocity.x() % state.next.velocity.y() % state.next.velocity.z() %
              state.next.pressure % state.next.energy;
        return os;
    }

    /// @brief Получить экземпляр расширенного вектора состояния
    [[nodiscard]] static State datatype();

    /// @brief Конструктор класса
    explicit SmFluid(const phys::Eos &eos, Fluxes flux = Fluxes::HLLC);

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

    void check_asserts(Mesh &mesh, std::string msg);

    /// @brief Установить флаги адаптации
    void set_flags(Mesh& mesh);

    /// @brief Распределитель данных при адаптации
    Distributor distributor() const;

    void set_accuracy(int acc);

    template<typename Test>
    void init_cells(Mesh &mesh, const Test &test) {
        static const SmFluid::State U = SmFluid::datatype();
        // Заполняем начальные данные
        for (auto cell: mesh) {
            cell(U).rho = test.density(cell.center());
            cell(U).v = test.velocity(cell.center());
            cell(U).p = test.pressure(cell.center());
            cell(U).e = m_eos.energy_rp(cell(U).rho, cell(U).p);
        }
    }

    double get_m_dt() {
        return m_dt;
    }

private:
    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    double compute_dt(Mesh &mesh);

    /// @brief Расчёт потоков
    void fluxes(Mesh &mesh);

    /// @brief Обновление ячеек
    void swap(Mesh &mesh);

    void compute_grad(Mesh &mesh,  const std::function<smf::PState(Cell &)> &to_state);

    void fluxes_stage1(Mesh &mesh);

    void fluxes_stage2(Mesh &mesh);

    smf::Flux calc_flux_extra(Cell &cell, bool from_begin);

public:

    ~SmFluid() = default;

protected:
    const phys::Eos &m_eos;
    NumFlux::Ptr m_nf; ///< Метод расчёта потока
    int m_acc = 1;
    double m_time = 0.0; ///< Прошедшее время
    size_t m_step = 0; ///< Количество шагов расчёта
    double m_CFL; ///< Число Куранта
    double m_dt; ///< Шаг интегрирования
};

} // namespace zephyr