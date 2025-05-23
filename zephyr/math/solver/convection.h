#pragma once

#include <zephyr/mesh/euler/soa_mesh.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math {

using zephyr::geom::Vector3d;
using zephyr::mesh::QCell;
using zephyr::mesh::SoaMesh;
using zephyr::mesh::Storable;
using zephyr::mesh::Distributor;

/// @brief Класс для моделирования уравнения переноса.
/// Пример образцово показательного решателя.
class Convection {
public:
    /// @brief Расширенный вектор состояния на котором решается задача
    Storable<double> u_curr, du_dx, du_dy, du_dz, u_half, u_next;

    /// @brief Конструктор класса, параметры по умолчанию:
    /// 3 порядок, CFL = 0.5, ограничитель minmod.
    Convection();

    /// @brief Виртуальный деструктор, наследование используется
    /// для переопределения функции velocity.
    virtual ~Convection() = default;

    /// @brief Добавить типы для хранения на сетку
    void add_types(SoaMesh& mesh);

    /// @brief Число Куранта
    double CFL() const;

    /// @brief Установить число Куранта
    void set_CFL(double C);

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double dt() const;

    /// @brief Название используемого ограничителя
    std::string limiter() const;

    /// @brief Установить ограничитель по названию
    void set_limiter(const std::string& limiter);

    /// @brief Используемый порядок точности
    int accuracy() const;

    /// @brief Установить порядок точности
    void set_accuracy(int acc);

    /// @brief Векторное поле скорости
    /// @details Виртуальная функция, следует унаследоваться от класса
    /// Convection и написать собственную функцию скорости
    virtual Vector3d velocity(const Vector3d& c) const;

    /// @brief Один шаг интегрирования по времени
    void update(SoaMesh& mesh);

    /// @brief Установить флаги адаптации
    void set_flags(SoaMesh& mesh);

    /// @brief Распределитель данных при адаптации
    Distributor distributor() const;

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    double compute_dt(QCell& cell) const;

    /// @brief Посчитать градиент хранимой функции
    /// @param stage При stage = 0 производные считаются для предыдущего
    /// временного слоя, при stage = 1 производные считаются для
    /// промежуточного временного слоя.
    void compute_grad(QCell& cell, int stage) const;

    /// @brief Просуммировать потоки.
    /// @param stage При stage = 0 выполняется шаг предиктора, при stage = 1
    /// выполняется шаг коррекции.
    void fluxes(QCell& cell, int stage);

protected:
    int m_accuracy;     ///< Точность схемы (1/2/3)
    double m_CFL;       ///< Число Куранта
    double m_dt;        ///< Шаг интегрирования
    Limiter m_limiter;  ///< Ограничитель
};

} // namespace zephyr::mesh