#pragma once

#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/mesh/euler/eu_mesh.h>

namespace zephyr::geom {

using mesh::EuCell;
using mesh::EuMesh;
using mesh::AmrCells;
using mesh::Storable;
using geom::Vector3d;

/// @brief Восстановление интерфейса для многоматериальных задач.
/// Работает на сетках с произвольными данными, единственное условие: 
/// в данных должны быть поля для хранения трех величин
/// 1. Объемная доля
/// 2. Нормаль интерфейса в ячейке
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
    void update(EuMesh& mesh, int smoothing) const;

    /// @brief Преобразовать реконструкцию в сетку
    EuMesh body(EuMesh& mesh) const;

private:

    /// @brief Оценка нормали по производным объемной доли
    void compute_normal(EuCell& cell) const;

    void compute_normals(EuMesh& mesh) const;

    /// @brief Найти сечение ячейки по существующим нормалям
    void find_section(EuCell& cell) const;

    void find_sections(EuMesh& mesh) const;

    /// @brief Поправить нормали с учетом соседей
    void adjust_normal(EuCell& cell) const;

    void adjust_normals(EuMesh& mesh) const;

private:
    Storable<double> a;    ///< Тип извлечения объемной доли
    Storable<Vector3d> n;  ///< Тип для извлечения нормали
    Storable<Vector3d> p;  ///< Тип для извлечения точки сечения
};

} // namespace zephyr::geom