#pragma once

#include <zephyr/mesh/euler/eu_mesh.h>

namespace zephyr { namespace math {

using zephyr::mesh::EuCell;
using zephyr::mesh::EuMesh;
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
        Vector3d n2;    ///< Нормаль (вспомогательное значение)

        inline bool zero_normal() const {
            return n.squaredNorm() == 0.0;
        }
    };

    /// @brief Получить экземпляр расширенного вектора состояния
    static State datatype();

    /// @brief Конструктор класса, по умолчанию CFL = 0.5
    Transfer();

    /// @brief Число Куранта
    double CFL() const;

    /// @brief Установить число Куранта
    void set_CFL(double C);

    /// @brief Версия функции update
    void set_version(int ver);

    /// @brief Использовать расщепление по направлениям
    void dir_splitting(bool flag);

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double dt() const;

    /// @brief Векторное поле скорости
    /// @details Виртуальная функция, следует унаследоваться от класса
    /// Transfer и написать собственную функцию скорости
    virtual Vector3d velocity(const Vector3d& c) const;

    /// @brief Один шаг интегрирования по времени
    void update(EuMesh& mesh);

    /// @brief Обновить флаг направления
    void update_dir();

    void compute_normals(EuMesh& mesh, int smoothing = 0);

    void find_sections(EuMesh& mesh);

    /// @brief Установить флаги адаптации
    void set_flags(EuMesh& mesh);

    /// @brief Распределитель данных при адаптации
    Distributor distributor() const;

    AmrStorage body(EuMesh& mesh);

    AmrStorage scheme(EuMesh& mesh);

protected:

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта
    double compute_dt(EuCell& cell);

    double compute_dt(EuMesh& mesh);

    void find_section(EuCell& cell);

    /// @brief Посчитать нормали
    void compute_normal(EuCell& cell);

    void update_ver1(EuMesh& mesh);

    void update_ver2(EuMesh& mesh);

    void update_ver3(EuMesh& mesh);

    /// @brief Потоки по схеме CRP
    void fluxes_CRP(EuCell& cell, Direction dir = Direction::ANY);

    /// @brief Потоки по аналогу VOF
    void fluxes_VOF(EuCell& cell, Direction dir = Direction::ANY);

    /// @brief Потоки по схеме CRP, но a_sig выбирается по аналогу VOF
    void fluxes_MIX(EuCell& cell, Direction dir = Direction::ANY);

protected:

    void smooth_normal(EuCell& cell);

    void update_normal(EuCell& cell);


    double m_dt;      ///< Шаг интегрирования
    double m_CFL;     ///< Число Куранта
    int    m_ver;     ///< Версия функции update
    Direction m_dir;  ///< Направление на текущем шаге
};

} // namespace math
} // namespace zephyr