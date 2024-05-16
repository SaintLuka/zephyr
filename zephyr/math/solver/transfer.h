#pragma once

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/geom/interface_recovery.h>

namespace zephyr::math {

using zephyr::mesh::EuCell;
using zephyr::mesh::EuMesh;
using zephyr::mesh::AmrStorage;
using zephyr::mesh::Distributor;
using zephyr::geom::Vector3d;
using zephyr::geom::Direction;
using zephyr::geom::InterfaceRecovery;


/// @brief Класс для моделирования уравнения переноса с аналогом CRP.
class Transfer {
public:

    /// @brief Расширенный вектор состояния на котором решается задача
    struct State {
        double u1, u2;  ///< Объемные доли
        double du_dx;   ///< Производная
        double du_dy;   ///< Производная
        Vector3d n;     ///< Внешняя нормаль поверхности
        Vector3d p;     ///< Базисная точка поверхности
    };

    enum class Method {
        // Методики CRP с эвристическими формулами, допускают
        // расщепление по направлениям и расчеты на полигональной сетке
        CRP_V3,      ///<
        CRP_V5,      ///<
        CRP_SE,      ///< Формула Серёжкина
        CRP_N1,      ///< Формула с учетом нормалей
        CRP_N2,      ///< Формула с учетом нормалей

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

        // Методики с WENO интерполяцией, допускают расщепление по
        // направлениям, не подходят для полигональной сетки
        WENO,        ///< WENO интеполяция на грань
        WENO_CRP,    ///< WENO с CRP ограничением

        // Другие, название?
        OTHER
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
    void set_method(Method method);

    /// @brief Использовать расщепление по направлениям
    void dir_splitting(bool flag);

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
    double compute_dt(EuMesh& mesh);

    /// @brief Один шаг интегрирования по времени
    void update(EuMesh& mesh, Direction dir = Direction::ANY);

    /// @brief Подсеточная реконструкция границы
    /// @param smoothing Число итераций сглаживания
    void update_interface(EuMesh& mesh, int smoothing = 3);

    /// @brief Установить флаги адаптации
    void set_flags(EuMesh& mesh);

    /// @brief Распределитель данных при адаптации
    Distributor distributor() const;

    AmrStorage body(EuMesh& mesh);

    AmrStorage scheme(EuMesh& mesh);

protected:

    /// @brief Посчитать шаг интегрирования по времени с учетом
    /// условия Куранта (для одной ячейки)
    double compute_dt(EuCell& cell);

    void compute_slopes(EuMesh& mesh);


    void update_CRP(EuMesh& mesh, Direction dir);

    void update_VOF(EuMesh& mesh, Direction dir);

    void update_MUSCL(EuMesh& mesh, Direction dir);

    void update_WENO(EuMesh& mesh, Direction dir);

    void update_other(EuMesh& mesh, Direction dir);


    /// @brief Потоки по схеме CRP
    void fluxes_CRP(EuCell& cell, Direction dir = Direction::ANY);

    /// @brief Потоки по аналогу VOF
    void fluxes_VOF(EuCell& cell, Direction dir = Direction::ANY);

    /// @brief Потоки по схеме MUSCL
    void fluxes_MUSCL(EuCell& cell, Direction dir = Direction::ANY);


protected:

    double m_dt;      ///< Шаг интегрирования
    double m_CFL;     ///< Число Куранта
    Method m_method;  ///< Методика вычисления потоков

    InterfaceRecovery interface; ///< Реконструкция границы
};

} // namespace zephyr::math