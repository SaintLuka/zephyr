#pragma once

#include <vector>

#include <zephyr/geom/vector.h>
#include <zephyr/mesh/index.h>
#include <zephyr/mesh/memory.h>
#include <zephyr/utils/range.h>

// forward declaration для классов из geom
namespace zephyr::geom {
class Quad;
class SqQuad;
class Cube;
class SqCube;
}

namespace zephyr::mesh {

/// @brief Квадратичное отображение на квадрат/куб в зависимости от размерности
template <int dim>
using SqMap = std::conditional_t<dim < 3, geom::SqQuad, geom::SqCube>;

/// @brief Набор дублирующихся вершин ячеек в форме Structure of Arrays (набор массивов).
///
class AmrVerts final {
    // aliases inside class
    using Vector3d = geom::Vector3d;

    /// @brief Используются уникальные узлы? Если unique = false, тогда массивы
    /// index и alien пустые. В обратном случае все массивы одного размера.
    bool m_unique = false;

public:
    /// @brief Индексы первых вершин ячеек (CSR-структура)
    std::vector<index_t>  offsets = {0};

    /// @brief Координаты вершин (с дубликатами)
    std::vector<Vector3d> coords;

    /// @brief Индекс узла в массиве local_nodes (в реальном локальном
    /// хранилище или в удаленном)
    std::vector<index_t> index;

    /// @brief Индекс узла в массиве ghost_nodes
    /// (или -1, если узел с данного процесса).
    std::vector<index_t> ghost;


    /// @brief Пустые массивы по умолчанию
    AmrVerts() = default;

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
    void init_unique(index_t idx, index_t gst);

    /// @brief Расход памяти
    memory_t memory_usage() const;

    /// @brief Число вершин, оно же максимальное, хранение неактуальных вершин
    /// не допускается.
    int count(index_t ic) const {
        return offsets[ic + 1] - offsets[ic];
    }

    /// @brief Число вершин, оно же максимальное, хранение неактуальных вершин
    /// сейчас не допускается.
    int max_count(index_t ic) const {
        return offsets[ic + 1] - offsets[ic];
    }

    /// @brief Полный диапазон вершин ячейки
    range_t<index_t> range(index_t ic) const {
        return std::views::iota(offsets[ic], offsets[ic + 1]);
    }

    /// @brief Указатель на первую вершину ячейки
    Vector3d* coords_data(index_t ic) {
        return coords.data() + offsets[ic];
    }

    /// @brief Константный указатель на первую вершину
    const Vector3d* coords_data(index_t ic) const {
        return coords.data() + offsets[ic];
    }

    /// @brief Ссылка на вершины в форме набора узлов квадратичного отображения
    template <int dim>
    SqMap<dim>& mapping(index_t ic) {
        return *reinterpret_cast<SqMap<dim>*>(coords_data(ic));
    }

    /// @brief Ссылка на вершины в форме набора узлов квадратичного отображения
    template <int dim>
    const SqMap<dim>& mapping(index_t ic) const {
        return *reinterpret_cast<const SqMap<dim>*>(coords_data(ic));
    }
};

} // namespace zephyr::mesh
