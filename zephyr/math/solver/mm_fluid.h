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
using zephyr::phys::Materials;

using namespace zephyr::math::mmf;

class MmFluid {
public:

    /// @brief Расширенный вектор состояния на котором решается задача
    struct State {
        double density;      ///< Плотность смеси
        Vector3d velocity;   ///< Равновесная скорость
        double pressure;     ///< Равновесное давление
        double energy;       ///< Равновесная энергия
        double temperature;  ///< Равновесная температура

        Fractions mass_frac; ///< Массовые доли веществ
        Fractions vol_frac;  ///< Объемные доли веществ

        PState half;         ///< Состояние на полушаге
        PState next;         ///< Состояние на следующем шаге

        /// @brief Градиент вектора состояния
        PState d_dx, d_dy, d_dz;


        bool is_bad1() const { return get_state().is_bad(); }

        bool is_bad2() const { return half.is_bad(); }

        bool is_bad3() const { return next.is_bad(); }

        bool is_bad() const { return is_bad1() || is_bad3(); }

        /// @brief Собрать вектор состояния на предыдущем шаге
        PState get_state() const {
            return PState(density, velocity, pressure, energy,
                          temperature, mass_frac, vol_frac);
        }

        /// @brief Установить вектор состояния на предыдущем шаге
        void set_state(const PState &z) {
            density     = z.density;
            velocity    = z.velocity;
            pressure    = z.pressure;
            energy      = z.energy;
            temperature = z.temperature;

            mass_frac   = z.mass_frac;
            vol_frac    = z.vol_frac;
        }
    };

    /// @brief Получить экземпляр расширенного вектора состояния
    static State datatype();

    /// @brief В поток вывода
    friend std::ostream &operator<<(std::ostream &os, const State &state);


    /// @brief Конструктор класса
    explicit MmFluid(const phys::Materials &eos);

    /// @brief Декструктор
    ~MmFluid() = default;

    /// @brief Установить число Куранта
    void set_CFL(double CFL);

    /// @brief Задать точность метода (1 или 2)
    void set_accuracy(int acc);

    /// @brief Установить метод
    void set_method(Fluxes method);

    /// @brief Задать ускорение свободного падения
    void set_gravity(double g);

    /// @brief Число Куранта
    double CFL() const;

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double dt() const;

    void update(Mesh &mesh);

    /// @brief Установить флаги адаптации
    void set_flags(Mesh &mesh);

    /// @brief Распределитель данных при адаптации
    Distributor distributor() const;

private:

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    void compute_dt(Mesh &mesh);

    /// @brief Расчёт потоков
    void fluxes(Mesh &mesh);

    /// @brief Обновление ячеек
    void swap(Mesh &mesh);

    void compute_grad(Mesh &mesh, const std::function<mmf::PState(Cell &)> &to_state);

    void fluxes_stage1(Mesh &mesh);

    void fluxes_stage2(Mesh &mesh);

    Flux calc_flux(const PState& zL, const PState& zR, double hL, double hR, double dt);

protected:
    Materials mixture;  ///< Список материалов (смесь)
    NumFlux::Ptr m_nf;  ///< Метод расчёта потока
    int m_acc = 1;      ///< Порядок точности
    double m_CFL;       ///< Число Куранта
    double m_g = 0.0;   ///< Ускорение свободного падения, направлено против oy
    double m_dt;        ///< Шаг интегрирования
    bool m_crp;         ///< Композитная задача Римана
};

std::ostream &operator<<(std::ostream &os, const MmFluid::State &state) {
    os << boost::format(
            "Main state: density: %1%, velocity: {%2%, %3%, %4%}, pressure: %5%, temperature: %6%, energy: %7%, mass_frac: %8%\n") %
          state.density % state.velocity.x() % state.velocity.y() % state.velocity.z() % state.pressure %
          state.temperature % state.energy % state.mass_frac;
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

} // namespace zephyr