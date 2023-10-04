#pragma once

#include <zephyr/mesh/mesh.h>
#include <zephyr/math/cfd/limiter.h>
#include "zephyr/phys/eos/eos.h"
#include "zephyr/phys/tests/classic_test.h"
#include "riemann.h"
#include <zephyr/math/cfd/fluxes.h>

namespace zephyr { namespace math {

using zephyr::mesh::ICell;
using zephyr::mesh::Mesh;
using zephyr::mesh::Distributor;
using zephyr::geom::Vector3d;

class MmSolver {
public:

    /// @brief Расширенный вектор состояния на котором решается задача
    struct State {
        double rho1, rho2; ///< плотность
        Vector3d v1, v2; ///< скорость
        double p1, p2; ///< давление
        double e1, e2; ///< энергия
    };

    /// @brief Получить экземпляр расширенного вектора состояния
    [[nodiscard]] static State datatype();

    /// @brief Конструктор класса
    explicit MmSolver(const phys::Eos &eos, Fluxes flux);

    /// @brief Число Куранта
    [[nodiscard]] double CFL() const;

    /// @brief Установить число Куранта
    void set_CFL(double CFL);

    /// @brief Установить метод расчёта потока
    void set_num_flux(Fluxes flux);

    [[nodiscard]] double get_time() const;

    [[nodiscard]] size_t getStep() const;

    [[nodiscard]] std::string get_flux_name() const;

    /// @brief Шаг интегрирования на предыдущем вызове update()
    [[nodiscard]] double dt() const;

    /// @brief Проинициализировать значения в ячейках согласно тесту
    void init_cells(Mesh &mesh, const phys::ClassicTest &test);

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    double compute_dt(Mesh &mesh);

    /// @brief Расчёт потоков
    void fluxes(Mesh &mesh);

    /// @brief Обновление ячеек
    void update(Mesh &mesh);

    ~MmSolver() = default;

protected:
    const phys::Eos &m_eos;
    NumFlux::Ptr m_nf; ///< Метод расчёта потока
    double m_time = 0.0; ///< Прошедшее время
    size_t m_step = 0; ///< Количество шагов расчёта
    double m_CFL; ///< Число Куранта
    double m_dt; ///< Шаг интегрирования
};

} // namespace math
} // namespace zephyr