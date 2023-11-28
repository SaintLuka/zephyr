/// @file Простое, но очень полезное определеине.

#include <zephyr/geom/primitives/mov_node.h>
#include <zephyr/geom/primitives/mov_cell.h>
#include <zephyr/mesh/storage.h>

namespace zephyr::mesh {

/// @brief Хранилище с подвижными ячейками
using CellStorage = Storage<zephyr::geom::MovCell>;

/// @brief Хранилище с узлами подвижной сетки
using NodeStorage = Storage<zephyr::geom::MovNode>;

} // namespace zephyr::mesh