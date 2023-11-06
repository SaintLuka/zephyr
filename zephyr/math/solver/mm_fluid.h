#pragma once

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
        double rho1, rho2; ///< плотность
        Vector3d v1, v2; ///< скорость
        double p1, p2; ///< давление
        double e1, e2; ///< энергия
        double t1, t2; ///< температура
        mmf::Fractions mass_frac1, mass_frac2; ///< доли веществ

        [[nodiscard]] bool is_bad2() const {
            return isinf(rho2) || isnan(rho2) ||
                   isinf(v2.x()) || isnan(v2.x()) ||
                   isinf(v2.y()) || isnan(v2.y()) ||
                   isinf(v2.z()) || isnan(v2.z()) ||
                   isinf(p2) || isnan(p2) ||
                   isinf(e2) || isnan(e2) ||
                   isinf(t2) || isnan(t2) ||
                   mass_frac2.empty();
        }

        [[nodiscard]] bool is_bad1() const {
            return isinf(rho1) || isnan(rho1) ||
                   isinf(v1.x()) || isnan(v1.x()) ||
                   isinf(v1.y()) || isnan(v1.y()) ||
                   isinf(v1.z()) || isnan(v1.z()) ||
                   isinf(p1) || isnan(p1) ||
                   isinf(e1) || isnan(e1) ||
                   isinf(t1) || isnan(t1) ||
                   mass_frac1.empty();
        }

        [[nodiscard]] bool is_bad() const {
            return is_bad1() || is_bad2();
        }

        [[nodiscard]] mmf::PState to_pstate() const {
            return mmf::PState(p1, v1, p1, e1, t1, mass_frac1);
        }
    };

    friend std::ostream &operator<<(std::ostream &os, const State &state) {
        os << boost::format(
                "State 1: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, temperature: %6%, energy: %7%, mass_frac: %8%\n") %
              state.rho1 % state.v1.x() % state.v1.y() % state.v1.z() % state.p1 %
              state.t1 % state.e1 % state.mass_frac1;
        os << boost::format(
                "State 2: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, temperature: %6%, energy: %7%, mass_frac: %8%\n") %
              state.rho2 % state.v2.x() % state.v2.y() % state.v2.z() % state.p2 %
              state.t2 % state.e2 % state.mass_frac2;
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

    /// @brief Установить метод расчёта потока
    void set_num_flux(Fluxes flux);

    [[nodiscard]] double get_time() const;

    [[nodiscard]] size_t get_step() const;

    [[nodiscard]] std::string get_flux_name() const;

    /// @brief Шаг интегрирования на предыдущем вызове update()
    [[nodiscard]] double dt() const;

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    double compute_dt(Mesh &mesh);

    /// @brief Расчёт потоков
    void fluxes(Mesh &mesh);

    /// @brief Обновление ячеек
    void update(Mesh &mesh);

    ~MmFluid() = default;

protected:
    const phys::Materials mixture;
    NumFlux::Ptr m_nf; ///< Метод расчёта потока
    double m_time = 0.0; ///< Прошедшее время
    size_t m_step = 0; ///< Количество шагов расчёта
    double m_CFL; ///< Число Куранта
    double m_dt; ///< Шаг интегрирования
};

} // namespace zephyr