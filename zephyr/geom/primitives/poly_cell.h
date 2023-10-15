#pragma once

#include <zephyr/geom/primitives/poly_faces.h>
#include <zephyr/geom/primitives/poly_nodes.h>
#include <zephyr/geom/primitives/element.h>

namespace zephyr::geom {

/// @struct Обязательные данные ячейки сетки
class PolyCell : public Element {
public:
    // Геометрия ячейки

    int       dim;     ///< Размерность ячейки
    Vector3d  coords;  ///< Барицентр ячейки
    double    size;    ///< Линейный размер ячейки
    PolyNodes nodes;   ///< Вершины ячейки
    PolyFaces faces;   ///< Список граней


};

} // namespace zephyr::geom