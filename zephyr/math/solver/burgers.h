#pragma once

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math {

using mesh::EuMesh;
using mesh::EuCell;
using mesh::Storable;
using mesh::Distributor;
using geom::Vector3d;

/// @brief Решение невязкого многомерного уравнения Бюргерса (уравение Хопфа)
/// d/dt (u_i) + u_k du_i / dx_k = 0
class BurgersSolver {
public:

    /// @brief Расширенный вектор состояния на котором решается задача
    struct Parts {
        Storable<Vector3d> init;  ///< Состояние на основном слое
        Storable<Vector3d> half;  ///< Состояние на полушаге
        Storable<Vector3d> next;  ///< Состояние на следующем шаге

        /// @brief Градиент вектора состояния
        Storable<Vector3d> d_dx, d_dy, d_dz;
    };

    Parts part;

    /// @brief Конструктор класса
    BurgersSolver();

    /// @brief Добавить типы на сетку
    Parts add_types(EuMesh& mesh);

    /// @brief Установить число Куранта
    void set_CFL(double CFL);

    /// @brief Задать точность метода (1 или 2)
    void set_accuracy(int acc);

    /// @brief Установить осевую симметрию
    void set_axial(bool axial = true);

    /// @brief Установить ограничитель градиента
    void set_limiter(const std::string& limiter);

    /// @brief Число Куранта
    double CFL() const;

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double dt() const;

    /// @brief Установить шаг интегрирования по времени
    void set_max_dt(double dt);

    /// @brief Выполнить шаг интегрирования по времени
    void update(EuMesh &mesh);

    /// @brief Установить флаги адаптации
    void set_flags(EuMesh& mesh) const;

    /// @brief Распределитель данных при адаптации
    /// @param type Тип "const" или "slope" переноса при разбиении
    Distributor distributor(const std::string& type = "slope") const;


    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    void compute_dt(EuMesh &mesh);

    /// @brief Расчёт потоков
    void fluxes(EuMesh &mesh) const;

    /// @brief Обновление ячеек
    void swap(EuMesh &mesh) const;

    /// @brief Вычислить производные
    void compute_grad(EuMesh &mesh) const;

    /// @brief Вычислить потоки на стадии предиктора
    void fluxes_stage1(EuMesh &mesh) const;

    /// @brief Вычислить потоки на стадии корректора
    void fluxes_stage2(EuMesh &mesh) const;

protected:
    int m_acc = 1;           ///< Порядок точности
    Limiter m_limiter;       ///< Ограничитель градиента
    double m_CFL = 0.5;      ///< Число Куранта
    double m_dt = NAN;       ///< Шаг интегрирования
    double m_max_dt=1.e300;  ///< Максимальный шаг интегрирования
};

} // namespace zephyr::math