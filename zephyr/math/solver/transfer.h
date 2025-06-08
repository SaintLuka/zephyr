#pragma once

#include <zephyr/mesh/euler/soa_mesh.h>
#include <zephyr/geom/interface_recovery.h>
#include <zephyr/math/cfd/limiter.h>

namespace zephyr::math {

using zephyr::geom::Vector3d;
using zephyr::geom::InterfaceRecovery;

using zephyr::mesh::QCell;
using zephyr::mesh::SoaMesh;
using zephyr::mesh::AmrCells;
using zephyr::mesh::Storable;
using zephyr::mesh::Direction;
using zephyr::mesh::AmrStorage;
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

        // Методики типа VOF с подсеточной реконструкицей границ,
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
        WENO,        ///< WENO интеполяция на грань
        WENO_CRP,    ///< WENO с CRP ограничением
    };

    // Расширенный вектор состояния на котором решается задача
    struct State {
        Storable<double> u1, u2;  ///< Объемные доли
        Storable<Vector3d> n;     ///< Внешняя нормаль поверхности
        Storable<Vector3d> p;     ///< Базисная точка поверхности

        // Градиенты нужны для схемы MUSCL
        Storable<double> du_dx;
        Storable<double> du_dy;
    };

    // Доступ к данным в хранилище
    State data;


    /// @brief Конструктор класса, по умолчанию CFL = 0.5
    Transfer();

    /// @brief Добавить типы для хранения на сетку
    void add_types(SoaMesh& mesh);

    /// @brief Число Куранта
    double CFL() const;

    /// @brief Установить число Куранта
    void set_CFL(double C);

    /// @brief Расчетный метод
    Method method() const;

    /// @brief Версия функции update
    void set_method(Method method);

    /// @brief Шаг интегрирования на предыдущем вызове update()
    double get_dt() const;

    /// @brief Установить временной шаг
    void set_dt(double dt);

    /// @brief Векторное поле скорости
    /// @details Виртуальная функция, следует унаследоваться от класса
    /// Transfer и написать собственную функцию скорости
    virtual Vector3d velocity(const Vector3d& c) const;

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта (для всех ячеек)
    double compute_dt(SoaMesh& mesh);

    /// @brief Один шаг интегрирования по времени
    void update(SoaMesh& mesh, Direction dir = Direction::ANY);

    /// @brief Подсеточная реконструкция границы
    /// @param smoothing Число итераций сглаживания
    void update_interface(SoaMesh& mesh, int smoothing = 3);

    /// @brief Установить флаги адаптации
    void set_flags(SoaMesh& mesh);

    /// @brief Распределитель данных при адаптации
    Distributor distributor() const;

    SoaMesh body(SoaMesh& mesh);

protected:

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта (для одной ячейки)
    double compute_dt(QCell& cell);

    void compute_slopes(SoaMesh& mesh) const;

    void update_CRP(SoaMesh& mesh, Direction dir);

    void update_VOF(SoaMesh& mesh, Direction dir);

    void update_MUSCL(SoaMesh& mesh, Direction dir);

    void update_WENO(SoaMesh& mesh, Direction dir);


    /// @brief Потоки по схеме CRP
    void fluxes_CRP(QCell& cell, Direction dir = Direction::ANY);

    /// @brief Потоки по аналогу VOF
    void fluxes_VOF(QCell& cell, Direction dir = Direction::ANY);

    /// @brief Потоки по схеме MUSCL
    void fluxes_MUSCL(QCell& cell, Direction dir = Direction::ANY);

protected:
    double m_dt;       ///< Шаг интегрирования
    double m_CFL;      ///< Число Куранта
    Method m_method;   ///< Методика вычисления потоков
    Limiter m_limiter; ///< Ограничитель для MUSCL_MC

    /// @brief Реконструкция границы
    InterfaceRecovery interface;
};

} // namespace zephyr::math