#pragma once

namespace zephyr::geom {

/// @brief Поддерживаемые типы ячеек сетки
/// Повторяют VTK типы с добавлением типов AMR2D, AMR3D, POLYHEDRON
enum class CellType : int {
    TRIANGLE,    ///< Треугольник
    QUAD,        ///< Четырехугольник
    POLYGON,     ///< Произвольный полигон
    AMR2D,       ///< Двумерная AMR-ячейка
    TETRA,       ///< Тетраэдр
    PYRAMID,     ///< Пирамида с четырехугольным основанием
    WEDGE,       ///< Призма с треугольным основанием
    HEXAHEDRON,  ///< Топологический куб
    POLYHEDRON,  ///< Произвольный многогранник
    AMR3D        ///< Трехмерная AMR-ячейка
};

} // namespace zephyr::geom