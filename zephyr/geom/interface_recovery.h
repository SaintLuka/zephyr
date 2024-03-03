#pragma once

#include <zephyr/mesh/euler/eu_mesh.h>

namespace zephyr::geom {

using zephyr::mesh::EuCell;
using zephyr::mesh::EuMesh;
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

    InterfaceRecovery(int volume_offset, int normal_offset, int origin_offset)
        : a{volume_offset}, n{normal_offset}, p{origin_offset} { }

    /// @brief Обновить подсеточную реконструкцию
    void update(EuMesh& mesh, int smoothing) const;

    /// @brief Преобразовать реконструкцию в сетку
    AmrStorage body(EuMesh& mesh) const;

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
    VarExtra<double> a;    ///< Тип извлечения объемной доли
    VarExtra<Vector3d> n;  ///< Тип для извлечения нормали
    VarExtra<Vector3d> p;  ///< Тип для извлечения точки сечения
};

} // namespace zephyr::geom