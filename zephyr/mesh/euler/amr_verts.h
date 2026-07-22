#pragma once

#include <vector>

#include <zephyr/geom/vector.h>
#include <zephyr/mesh/index.h>
#include <zephyr/mesh/memory.h>


namespace zephyr::mesh {

/// @brief Набор дублирующихся вершин ячеек в форме Structure of Arrays (набор массивов).
///
class AmrVerts final {
    // aliases inside class
    using Vector3d = geom::Vector3d;

    /// @brief Используются уникальные узлы? Если unique = false, тогда массивы
    /// index и alien пустые. В обратном случае все массивы одного размера.
    bool m_unique = false;

public:
    /// @brief Координаты вершин (с дубликатами)
    std::vector<Vector3d> coords;

    /// @brief Индекс узла в массиве local_nodes (в реальном локальном
    /// хранилище или в удаленном)
    std::vector<index_t> index;

    /// @brief Индекс узла в массиве ghost_nodes
    /// (или -1, если узел с данного процесса).
    std::vector<index_t> ghost;


    /// @brief Используются уникальные узлы?
    bool unique() const { return m_unique; }

    /// @brief Число вершин
    index_t size() const { return coords.size(); }

    /// @brief Изменить размер под число вершин
    void resize(index_t n_verts);

    /// @brief Расширить буфер под число вершин
    void reserve(index_t n_verts);

    /// @brief Сжать до актуальных размеров
    void shrink_to_fit();

    /// @brief Оператор доступа к координате
    Vector3d& operator[](index_t idx) { return coords[idx]; }

    /// @brief Оператор доступа к координате
    const Vector3d& operator[](index_t idx) const { return coords[idx]; }

    /// @brief Забыть об уникальных узлах
    void clear_unique();

    /// @brief Инициализировать массивы для уникальных узлов
    void init_unique();

    /// @brief Расход памяти
    memory_t memory_usage() const;
};

} // namespace zephyr::mesh
