/// @file mov_storage.h
/// @brief Простое, но очень полезное определение

#include <zephyr/mesh/primitives/mov_node.h>
#include <zephyr/mesh/primitives/mov_cell.h>
#include <zephyr/mesh/storage.h>

namespace zephyr::mesh {

/// @brief Хранилище с подвижными ячейками
using CellStorage = Storage<MovCell>;

/// @brief Хранилище с узлами подвижной сетки
using NodeStorage = Storage<MovNode>;

} // namespace zephyr::mesh