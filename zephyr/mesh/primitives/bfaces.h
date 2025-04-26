#pragma once

#include <zephyr/geom/cell_type.h>
#include <zephyr/mesh/primitives/Side3D.h>
#include <zephyr/mesh/primitives/bface.h>

namespace zephyr::mesh {

/// @brief Массив граней ячейки
class BFaces {
    using CellType = zephyr::geom::CellType;

public:
    /// @brief Максимальное число граней
    static const int max_count = 24;

    /// @brief Конструктор инициализирует индексы граней
    /// @param count Число вершин, требуется при инициализации полигонов
    explicit BFaces(CellType ctype, int count = -1);

    /// @brief Число актуальных граней в списке
    int count() const;

    /// @brief Доступ к грани по индексу
    inline BFace &operator[](int idx) {
        return m_faces[idx];
    }

    /// @brief Доступ к грани по индексу
    inline const BFace &operator[](int idx) const {
        return m_faces[idx];
    }

    /// @brief Начало списка
    inline std::array<BFace, max_count>::iterator begin() {
        return m_faces.begin();
    }

    /// @brief Начало списка
    inline std::array<BFace, max_count>::const_iterator begin() const {
        return m_faces.begin();
    }

    /// @brief Конец списка
    inline std::array<BFace, max_count>::iterator end() {
        return m_faces.end();
    }

    /// @brief Конец списка
    inline std::array<BFace, max_count>::const_iterator end() const {
        return m_faces.end();
    }

private:
    /// @brief Массив граней ячейки
    std::array<BFace, max_count> m_faces;
};

} // namespace zephyr::mesh