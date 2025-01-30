#pragma once

#include <zephyr/geom/vector.h>
#include <zephyr/geom/generator/rectangle.h>

#include <zephyr/mesh/euler/eu_cell.h>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/fractions.h>
#include <zephyr/phys/matter/materials.h>

namespace zephyr::phys {

using zephyr::geom::Vector3d;
using zephyr::geom::Boundary;
using zephyr::geom::generator::Rectangle;

/// @class Абстрактный класс для задания начальных данных.
/// Самое общее задание начальных данных из всех возможных.
/// для некоторых тестов доступно точное решение.
class IVP {
protected:
    Materials m_materials;  ///< Список материалов

public:
    /// @brief Получить название теста
    virtual std::string name() const { return "SomeTest"; }

    /// @brief Конечный момент времени
    virtual double max_time() const = 0;

    /// @brief Список материалов
    const Materials& materials() const;

    /// @brief Смесь PT-замыкание
    MixturePT mixture_PT() const;

    /// @brief Получить уравнение состояния по индексу
    Eos::Ptr get_eos(int idx = 0) const;


    // Распределения основных параметров в начальный момент времени
    // в зависимости от координаты

    /// @brief Индекс материала в точке
    virtual int index(const Vector3d &r) const { return 0; }

    /// @brief Начальная плотность
    virtual double density(const Vector3d &r) const = 0;

    /// @brief Начальная скорость
    virtual Vector3d velocity(const Vector3d &r) const = 0;

    /// @brief Начальное давление
    virtual double pressure(const Vector3d &r) const = 0;

    /// @brief Начальная внутренняя энергия (по умолчанию УрС)
    virtual double energy(const Vector3d &r) const;

    /// @brief Начальная температура (по умолчанию УрС)
    virtual double temperature(const Vector3d &r) const;

    /// @brief Характеристические функции компонент, "объемные доли"
    /// по умолчанию для одноматериальной
    virtual Fractions fractions(const Vector3d &r) const;

    /// @brief Характеристическая функция для материала с индексом idx
    /// по умолчанию true для idx = 0.
    bool inside(const Vector3d &r, int idx) const;

    /// @defgroup
    /// Эволюция параметров от времени. Функции следует определять
    /// для задач с точным решением. По умолчанию возращают NAN и нули.
    /// @{

    virtual int index_t(const Vector3d& r, double t) const;

    virtual double density_t(const Vector3d &r, double t) const;

    virtual Vector3d velocity_t(const Vector3d &r, double t) const;

    virtual double pressure_t(const Vector3d &r, double t) const;

    virtual double energy_t(const Vector3d &r, double t) const;

    virtual double temperature_t(const Vector3d& r, double t) const;

    virtual Fractions fractions_t(const Vector3d &r, double t) const;

    bool inside_t(const Vector3d &r, int idx, double t) const;

    /// @}


    /// @defgroup
    /// Подсеточное задание начальных параметров, при значении параметра
    /// n < 2 возвращает простое значение в центре ячейки.
    /// при n > 1 разбивает ячейку на n по каждой стороне и суммирует.
    /// @{

    /// @brief Плотность вещества ячейки
    double density_mean(mesh::EuCell& cell, int n) const;

    /// @brief Плотность момента импульса ячейки
    Vector3d momentum_mean(mesh::EuCell& cell, int n) const;

    /// @brief Удельная полная энергия ячейки
    double energy_mean(mesh::EuCell& cell, int n) const;

    /// @brief Массовые доли компонент смеси
    Fractions mass_fractions(mesh::EuCell& cell, int n) const;

    /// @brief Объемные доли компонент смеси
    Fractions volume_fractions(mesh::EuCell& cell, int n) const;

    /// @}
};

} // namespace zephyr::phys