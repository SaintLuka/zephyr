#pragma once

#include <zephyr/geom/polygon.h>

#include <zephyr/mesh/primitives/element.h>
#include <zephyr/mesh/primitives/bnodes.h>
#include <zephyr/mesh/primitives/bfaces.h>

namespace zephyr::mesh {

/// @struct Обязательные данные ячейки сетки
class MovCell : public Element {
    using Vector3d = zephyr::geom::Vector3d;
    using Polygon  = zephyr::geom::Polygon;
public:
    // Геометрия ячейки

    int      dim;     ///< Размерность ячейки
    Vector3d center;  ///< Барицентр ячейки
    double   size;    ///< Линейный размер ячейки
    BNodes   nodes;   ///< Вершины ячейки
    BFaces   faces;   ///< Список граней

    /// @brief Двумерная полигональная ячейка.
    /// @details Массив nodes не инициализируется
    explicit MovCell(const Polygon& poly);

    /// @brief Площадь (в 2D) или объем (в 3D) ячейки
    double volume() const;
};

} // namespace zephyr::mesh