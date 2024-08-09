#pragma once

#include <zephyr/mesh/mesh.h>
#include <zephyr/phys/eos/eos.h>
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
        double density;     ///< Плотность
        Vector3d velocity;  ///< Скорость
        double pressure;    ///< Давление
        double energy;      ///< Энергия

        PState half;        ///< Состояние на полушаге
        PState next;        ///< Состояние на следующем шаге

        /// @brief Градиент вектора состояния
        PState d_dx, d_dy, d_dz;

        /// @brief Собрать вектор состояния на предыдущем шаге
        PState get_state() const {
            return PState(density, velocity, pressure, energy);
        }

        /// @brief Установить вектор состояния на предыдущем шаге
        void set_state(const PState &z) {
            density  = z.density;
            velocity = z.velocity;
            pressure = z.pressure;
            energy   = z.energy;
        }
    };

    /// @brief Получить экземпляр расширенного вектора состояния
    static State datatype();

    /// @brief В поток вывода
    friend std::ostream &operator<<(std::ostream &os, const State &state);


    /// @brief Конструктор класса
    explicit SmFluid(const phys::Eos &eos);

    /// @brief Декструктор
    ~SmFluid() = default;

    /// @brief Установить число Куранта
    void set_CFL(double CFL);

    /// @brief Задать точность метода (1 или 2)
    void set_accuracy(int acc);

    /// @brief Установить метод
    void set_method(Fluxes method);

    /// @brief Число Куранта
    double CFL() const;

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double dt() const;

    /// @brief Выполнить шаг интегрирования по времени
    void update(Mesh &mesh);

    /// @brief Установить флаги адаптации
    void set_flags(Mesh& mesh);

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

    /// @brief Вычислить производные
    void compute_grad(Mesh &mesh,  const std::function<smf::PState(Cell &)> &get_state);

    /// @brief Вычислить потоки на стадии предиктора
    void fluxes_stage1(Mesh &mesh);

    /// @brief Вычислить потоки на стадии корректора
    void fluxes_stage2(Mesh &mesh);

protected:
    const phys::Eos &m_eos;  ///< Уравнение состояния
    NumFlux::Ptr m_nf;       ///< Метод расчёта потока
    int m_acc = 1;           ///< Порядок точности
    double m_CFL;            ///< Число Куранта
    double m_dt;             ///< Шаг интегрирования
};

std::ostream &operator<<(std::ostream &os, const SmFluid::State &state) {
    os << "State1: " << state.get_state() << "\n";
    os << "State2: " << state.half << "\n";
    return os;
}

} // namespace zephyr::math