#pragma once
#include <stdexcept>

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

inline int get_dimension(CellType type) {
    if (type == CellType::TRIANGLE ||
        type == CellType::QUAD ||
        type == CellType::POLYGON ||
        type == CellType::AMR2D) {
        return 2;
    }
    if (type == CellType::TETRA  ||
        type == CellType::PYRAMID ||
        type == CellType::WEDGE ||
        type == CellType::HEXAHEDRON ||
        type == CellType::POLYHEDRON ||
        type == CellType::AMR3D) {
        return 3;
    }
    throw std::runtime_error("dimension(CellType) error: unknown cell type");
}

} // namespace zephyr::geom