#pragma once

#include <zephyr/mesh/mesh.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr { namespace math {

using zephyr::mesh::ICell;
using zephyr::mesh::Mesh;
using zephyr::mesh::Distributor;
using zephyr::geom::Vector3d;

/// @brief Класс для моделирования уравнения переноса.
/// Пример образцово показательного решателя.
class Convection {
public:

    /// @brief Расширенный вектор состояния на котором решается задача
    struct State {
        double u1, ux, uy, uh, u2;
    };

    /// @brief Получить экземпляр расширенного вектора состояния
    static State datatype();

    /// @brief Конструктор класса, параметры по умолчанию:
    /// 3 порядок, CFL = 0.5, ограничитель minmod.
    Convection();

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
    void update(Mesh& mesh);

    /// @brief Установить флаги адаптации
    void set_flags(Mesh& mesh);

    /// @brief Распределитель данных при адаптации
    Distributor distributor() const;

protected:

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    double compute_dt(const ICell& cell) const;

    /// @brief Посчитать граниент хранимой функции
    /// @param stage При stage = 0 производные считаются для предыдущего
    /// временного слоя, при stage = 1 производные считаются для
    /// промежуточного временного слоя.
    void compute_grad(ICell& cell, int stage);

    /// @brief Просуммировать потоки.
    /// @param stage При stage = 0 выполняется шаг предиктора, при stage = 1
    /// выполняется шаг коррекции.
    void fluxes(ICell& cell, int stage);

protected:
    int m_accuracy;     ///< Точность схемы (1/2/3)
    double m_CFL;       ///< Число Куранта
    double m_dt;        ///< Шаг интегрирования
    Limiter m_limiter;  ///< Ограничитель
};

} // namespace math
} // namespace zephyr