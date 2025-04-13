#pragma once

#include <zephyr/mesh/mesh.h>
#include <zephyr/phys/matter/eos/eos.h>
#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math {

using zephyr::mesh::Cell;
using zephyr::mesh::Mesh;
using zephyr::mesh::Distributor;
using zephyr::geom::Vector3d;
using zephyr::phys::Eos;

using namespace smf;

/// @class SmFluid sm_fluid.h
/// @brief Single-Material Fluid. Решатель для классической
/// одноматериальной гидро- и газодинамики.
class SmFluidSoA {
public:

    struct State { };

    std::vector<PState> curr;
    std::vector<PState> half;
    std::vector<PState> next;
    std::vector<PState> d_dx;
    std::vector<PState> d_dy;
    std::vector<PState> d_dz;

    /// @brief Получить экземпляр расширенного вектора состояния
    static State datatype();

    /// @brief В поток вывода
    friend std::ostream &operator<<(std::ostream &os, const State &state);


    /// @brief Конструктор класса
    explicit SmFluidSoA(Eos::Ptr eos);

    /// @brief Декструктор
    ~SmFluidSoA() = default;

    /// @brief Установить число Куранта
    void set_CFL(double CFL);

    /// @brief Задать точность метода (1 или 2)
    void set_accuracy(int acc);

    /// @brief Установить осевую симметрию
    void set_axial(bool axial = true);

    /// @brief Установить метод
    void set_method(Fluxes method);

    /// @brief Установить ограничитель градиента
    void set_limiter(const std::string& limiter);

    /// @brief Число Куранта
    double CFL() const;

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double dt() const;

    /// @brief Установить шаг интегрирования по времени
    void set_max_dt(double dt);

    /// @brief Выполнить шаг интегрирования по времени
    void update(Mesh &mesh);

    /// @brief Установить флаги адаптации
    void set_flags(Mesh& mesh);

    /// @brief Распределитель данных при адаптации
    /// @param type Тип "const" или "slope" переноса при разбиении
    Distributor distributor(const std::string& type = "slope") const;


    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    void compute_dt(Mesh &mesh);

    /// @brief Расчёт потоков
    void fluxes(Mesh &mesh);

    /// @brief Обновление ячеек
    void swap(Mesh &mesh);

    /// @brief Вычислить производные
    void compute_grad(Mesh &mesh);

    /// @brief Вычислить производные
    void compute_grad(Mesh &mesh, const std::function<smf::PState(Cell &)> &get_state);

    /// @brief Вычислить потоки на стадии предиктора
    void fluxes_stage1(Mesh &mesh);

    /// @brief Вычислить потоки на стадии корректора
    void fluxes_stage2(Mesh &mesh);

protected:
    Eos::Ptr m_eos;          ///< Уравнение состояния
    NumFlux::Ptr m_nf;       ///< Метод расчёта потока
    int m_acc = 1;           ///< Порядок точности
    bool m_axial;            ///< Осевая симметрия
    Limiter m_limiter;       ///< Ограничитель градиента
    double m_CFL = 0.5;      ///< Число Куранта
    double m_dt;             ///< Шаг интегрирования
    double m_max_dt=1.e300;  ///< Максимальный шаг интегрирования
};

} // namespace zephyr::math