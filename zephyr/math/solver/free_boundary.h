#pragma once

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/phys/matter/eos/eos.h>
#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math {

using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::mesh::Storable;
using zephyr::geom::Vector3d;
using zephyr::phys::Eos;

using namespace smf;

/// @class FreeBoundary free_boundary.h
/// @brief Single-Material Fluid. Решатель для классической
/// одноматериальной гидро- и газодинамики.
class FreeBoundary {
public:

    /// @brief Расширенный вектор состояния на котором решается задача
    struct Parts {
        Storable<PState> init;  ///< Состояние на основном слое
        Storable<PState> next;  ///< Состояние на следующем шаге
    };

    Parts part;

    /// @brief Конструктор класса
    explicit FreeBoundary(Eos::Ptr eos);

    /// @brief Деструктор
    ~FreeBoundary() = default;

    /// @brief Добавить типы на сетку
    Parts add_types(EuMesh& mesh);

    /// @brief Установить число Куранта
    void set_CFL(double CFL);

    /// @brief Установить осевую симметрию
    void set_axial(bool axial = true);

    /// @brief Установить метод
    void set_method(Fluxes method);

    /// @brief Число Куранта
    double CFL() const;

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double dt() const;

    /// @brief Установить шаг интегрирования по времени
    void set_max_dt(double dt);

    /// @brief Выполнить шаг интегрирования по времени
    void update(EuMesh &mesh);

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    void compute_dt(EuMesh &mesh);

    /// @brief Расчёт потоков
    void fluxes(EuMesh &mesh) const;

    /// @brief Обновление ячеек
    void swap(EuMesh &mesh) const;

    /// @brief Вычислить производные
    void compute_grad(EuMesh &mesh) const;

protected:
    Eos::Ptr m_eos;          ///< Уравнение состояния
    int m_acc = 1;           ///< Порядок точности
    bool m_axial;            ///< Осевая симметрия
    double m_CFL = 0.5;      ///< Число Куранта
    double m_dt;             ///< Шаг интегрирования
    double m_max_dt=1.e300;  ///< Максимальный шаг интегрирования
};

} // namespace zephyr::math