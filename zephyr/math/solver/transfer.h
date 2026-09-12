#pragma once

#include <zephyr/mesh/mesh.h>
#include <zephyr/geom/plic.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math {

using zephyr::geom::Vector3d;
using zephyr::geom::Plic;

using zephyr::mesh::Cell;
using zephyr::mesh::Mesh;
using zephyr::mesh::RawCells;
using zephyr::mesh::Storable;
using zephyr::mesh::Direction;
using zephyr::mesh::Distributor;


/// @brief Класс для моделирования уравнения переноса с аналогом CRP.
class Transfer {
public:

    /// @brief Список методов решения
    enum class Method {
        // Методики CRP с эвристическими формулами, допускают
        // расщепление по направлениям и расчеты на полигональной сетке
        CRP_V3,      ///<
        CRP_V5,      ///<
        CRP_SE,      ///< Формула Серёжкина
        CRP_N1,      ///< Формула с учетом нормалей (точное пересечение)
        CRP_N2,      ///< Формула с учетом нормалей (average flux)

        // Методики типа VOF с подсеточной реконструкцией границ,
        // допускают расщепление по направлениям и расчеты
        // на полигональной сетке
        VOF,         ///< Обычный VOF
        VOF_CRP,     ///< VOF с CRP ограничением

        // Методики типа MUSCL, допускают расщепление по направлениям
        // и расчеты на полигональной сетке
        MUSCLd,      ///< MUSCL с расчетом производных
        MUSCLn,      ///< MUSCL с подсеточной реконструкцией
        MUSCLd_CRP,  ///< MUSCLd с CRP ограничением
        MUSCLn_CRP,  ///< MUSCLn с CRP ограничением
        MUSCL_MC,    ///< MUSCL с лимитированым градиентом (MC)
        MUSCL_MC_CRP,///< MUSCL_MC с CRP ограничением

        // Методики с WENO интерполяцией, допускают расщепление по
        // направлениям, не подходят для полигональной сетки
        WENO,        ///< WENO интерполяция на грань
        WENO_CRP,    ///< WENO с CRP ограничением
    };

    // Расширенный вектор состояния на котором решается задача
    struct State {
        Storable<double> u1, u2;  ///< Объемные доли
        Storable<double> p;       ///< Расстояние от центра ячейки до плоскости
        Storable<Vector3d> n;     ///< Внешняя нормаль поверхности

        // Градиент нужен для схемы MUSCL
        Storable<Vector3d> grad;
    };

    // Доступ к данным в хранилище
    State data;


    /// @brief Конструктор класса, по умолчанию CFL = 0.5
    Transfer();

    virtual ~Transfer() = default;

    /// @brief Задать размерность
    void set_dim(int dim);

    /// @brief Добавить типы для хранения на сетку
    State add_types(Mesh& mesh);

    /// @brief Число Куранта
    double CFL() const;

    /// @brief Установить число Куранта
    void set_CFL(double C);

    /// @brief Расчетный метод
    Method method() const;

    /// @brief Версия функции update
    void set_method(Method method);

    /// @brief Установить тип PLIC реконструкции
    void set_plic_type(Plic::Type type);

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double get_dt() const;

    /// @brief Установить временной шаг
    void set_dt(double dt);

    /// @brief Установить шаг интегрирования по времени
    void set_max_dt(double dt);

    /// @brief Векторное поле скорости
    /// @details Виртуальная функция, следует унаследоваться от класса
    /// Transfer и написать собственную функцию скорости
    virtual Vector3d velocity(const Vector3d& c) const;

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта (для всех ячеек)
    double compute_dt(Mesh& mesh) const;

    /// @brief Один шаг интегрирования по времени
    void update(Mesh& mesh, Direction dir = Direction::ANY);

    /// @brief Подсеточная реконструкция границы
    /// @param smoothing Число итераций сглаживания
    void update_interface(Mesh& mesh, int smoothing = 3) const;

    /// @brief Установить флаги адаптации
    void set_flags(Mesh& mesh) const;

    /// @brief Распределитель данных при адаптации
    Distributor distributor() const;

    Mesh body(Mesh& mesh) const;

protected:

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта (для одной ячейки)
    double compute_dt(Cell& cell) const;

    void compute_slopes(Mesh& mesh) const;

    void update_CRP(Mesh& mesh, Direction dir) const;

    void update_VOF(Mesh& mesh, Direction dir);

    void update_MUSCL(Mesh& mesh, Direction dir) const;

    void update_WENO(Mesh& mesh, Direction dir) const;


    /// @brief Потоки по схеме CRP
    void fluxes_CRP(Cell& cell, Direction dir = Direction::ANY) const;

    /// @brief Потоки по аналогу VOF
    void fluxes_VOF(Cell& cell, Direction dir = Direction::ANY) const;

    /// @brief Потоки по схеме MUSCL
    void fluxes_MUSCL(Cell& cell, Direction dir = Direction::ANY) const;

protected:
    double m_dt;       ///< Шаг интегрирования
    double m_CFL;      ///< Число Куранта
    int    m_dim;      ///< Размерность сетки/решателя
    Method m_method;   ///< Методика вычисления потоков
    Plic   m_plic;     ///< Метод реконструкции границы
    Limiter m_limiter; ///< Ограничитель для MUSCL_MC

    double m_max_dt=1.e300;  ///< Максимальный шаг интегрирования
};

} // namespace zephyr::math