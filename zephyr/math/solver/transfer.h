#pragma once

#include <zephyr/mesh/mesh.h>

namespace zephyr { namespace math {

using zephyr::mesh::ICell;
using zephyr::mesh::Mesh;
using zephyr::mesh::AmrStorage;
using zephyr::mesh::Distributor;
using zephyr::geom::Vector3d;
using zephyr::geom::Direction;

/// @brief Класс для моделирования уравнения переноса с аналогом CRP.
class Transfer {
public:

    /// @brief Расширенный вектор состояния на котором решается задача
    struct State {
        double u1, u2;  ///< Объемные доли
        Vector3d n;     ///< Внешняя нормаль поверхности
        Vector3d p;     ///< Базисная точка поверхности
    };

    /// @brief Получить экземпляр расширенного вектора состояния
    static State datatype();

    /// @brief Конструктор класса, по умолчанию CFL = 0.5
    Transfer();

    /// @brief Число Куранта
    double CFL() const;

    /// @brief Установить число Куранта
    void set_CFL(double C);

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double dt() const;

    /// @brief Векторное поле скорости
    /// @details Виртуальная функция, следует унаследоваться от класса
    /// Transfer и написать собственную функцию скорости
    virtual Vector3d velocity(const Vector3d& c) const;

    /// @brief Один шаг интегрирования по времени
    void update(Mesh& mesh, int ver = 1);

    void compute_normals(Mesh& mesh, int smoothing = 0);

    void find_sections(Mesh& mesh);

    /// @brief Установить флаги адаптации
    void set_flags(Mesh& mesh);

    /// @brief Распределитель данных при адаптации
    Distributor distributor() const;

    AmrStorage body(Mesh& mesh);

protected:

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    double compute_dt(ICell& cell);

    double compute_dt(Mesh& mesh);

    void find_section(ICell& cell);

    /// @brief Посчитать нормали
    void compute_normal(ICell& cell);

    void smooth_normal(ICell& cell);

    void update_ver1(Mesh& mesh);

    void update_ver2(Mesh& mesh);

    void update_ver3(Mesh& mesh);

    /// @brief Потоки по схеме CRP
    void fluxes_CRP(ICell& cell, Direction dir = Direction::ANY);

    /// @brief Потоки по аналогу VOF
    void fluxes_VOF(ICell& cell, Direction dir = Direction::ANY);

    /// @brief Потоки по схеме CRP, но a_sig выбирается по аналогу VOF
    void fluxes_MIX(ICell& cell, Direction dir = Direction::ANY);

protected:
    double m_CFL;       ///< Число Куранта
    double m_dt;        ///< Шаг интегрирования
};

} // namespace math
} // namespace zephyr