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

using namespace smf;

/// @class Single-Material Fluid.
/// @brief Класс решатель классической одноматериальной газодинамики
class SmFluid {
public:

    /// @brief Расширенный вектор состояния на котором решается задача
    struct State {
        double density;     ///< плотность
        Vector3d velocity;  ///< скорость
        double pressure;    ///< давление
        double energy;      ///< энергия

        PState half;        ///< Состояние на полушаге
        PState next;        ///< Состояние на следующем шаге

        /// @brief Градиент вектора состояния
        PState d_dx, d_dy, d_dz;

        /// @brief Собрать вектор состояния на предыдущем шаге
        smf::PState get_state() const {
            return smf::PState(density, velocity, pressure, energy);
        }

        /// @brief Установить вектор состояния на предыдущем шаге
        void set_state(const smf::PState &pstate) {
            density  = pstate.density;
            velocity = pstate.velocity;
            pressure = pstate.pressure;
            energy   = pstate.energy;
        }
    };

    /// @brief Получить экземпляр расширенного вектора состояния
    static State datatype();

    /// @brief В поток вывода
    friend std::ostream &operator<<(std::ostream &os, const State &state);


    /// @brief Конструктор класса
    explicit SmFluid(const phys::Eos &eos, Fluxes flux = Fluxes::HLLC);

    /// @brief Декструктор
    ~SmFluid() = default;

    /// @brief Число Куранта
    double CFL() const;

    /// @brief Установить число Куранта
    void set_CFL(double CFL);

    /// @brief Установить порядок точности
    void set_acc(int acc);

    /// @brief Установить метод
    void set_method(Fluxes method);

    double get_time() const;

    size_t get_step() const;

    std::string get_flux_name() const;

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double dt() const;

    void update(Mesh &mesh);

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
            cell(U).density  = test.density(cell.center());
            cell(U).velocity = test.velocity(cell.center());
            cell(U).pressure = test.pressure(cell.center());
            cell(U).energy   = m_eos.energy_rp(cell(U).density, cell(U).pressure);
        }
    }

private:
    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    void compute_dt(Mesh &mesh);

    /// @brief Расчёт потоков
    void fluxes(Mesh &mesh);

    /// @brief Обновление ячеек
    void swap(Mesh &mesh);

    void compute_grad(Mesh &mesh,  const std::function<smf::PState(Cell &)> &get_state);

    void fluxes_stage1(Mesh &mesh);

    void fluxes_stage2(Mesh &mesh);

protected:
    const phys::Eos &m_eos;  ///< Уравнение состояния
    NumFlux::Ptr m_nf;       ///< Метод расчёта потока
    int m_acc = 1;           ///< Порядок точности
    double m_CFL;            ///< Число Куранта
    double m_dt;             ///< Шаг интегрирования

    double m_time = 0.0;     ///< Прошедшее время
    size_t m_step = 0;       ///< Количество шагов расчёта
};

std::ostream &operator<<(std::ostream &os, const SmFluid::State &state) {
    os << "State1: " << state.get_state() << "\n";
    os << "State2: " << state.half << "\n";
    return os;
}

} // namespace zephyr::math