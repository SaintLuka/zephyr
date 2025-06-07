#pragma once

#include <zephyr/mesh/euler/soa_mesh.h>

namespace zephyr::geom {

using zephyr::mesh::QCell;
using zephyr::mesh::SoaMesh;
using zephyr::mesh::AmrCells;
using zephyr::mesh::Storable;
using zephyr::geom::Vector3d;
using zephyr::mesh::AmrStorage;
using zephyr::mesh::VarExtra;
using zephyr::mesh::Byte;


/// @brief Восстановление интерфейса для многоматериальных задач.
/// Работает на сетках с произвольными данными, единственное условие: 
/// в данных должны быть поля для хранения трех величин
/// 1. Объемная доля
/// 2. Номаль интерфейса в ячейке
/// 3. Точка плоскости внутри ячейке.
/// Последние две величины определяются с использованием данного класса.
class InterfaceRecovery {
public:

    InterfaceRecovery() = default;

    InterfaceRecovery(
            Storable<double> volume,
            Storable<Vector3d> normal,
            Storable<Vector3d> origin)
        : a{volume}, n{normal}, p{origin} { }

    /// @brief Обновить подсеточную реконструкцию
    void update(SoaMesh& mesh, int smoothing) const;

    /// @brief Преобразовать реконструкцию в сетку
    SoaMesh body(SoaMesh& mesh) const;

private:

    /// @brief Оценка нормали по производным объемной доли
    void compute_normal(QCell& cell) const;

    void compute_normals(SoaMesh& mesh) const;

    /// @brief Найти сечение ячейки по существующим нормалям
    void find_section(QCell& cell) const;

    void find_sections(SoaMesh& mesh) const;

    /// @brief Поправить нормали с учетом соседей
    void adjust_normal(QCell& cell) const;

    void adjust_normals(SoaMesh& mesh) const;

private:
    Storable<double> a;    ///< Тип извлечения объемной доли
    Storable<Vector3d> n;  ///< Тип для извлечения нормали
    Storable<Vector3d> p;  ///< Тип для извлечения точки сечения
};

} // namespace zephyr::geom