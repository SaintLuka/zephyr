/// @file Простое, но очень полезное определеине.

#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/mesh/storage.h>

namespace zephyr::mesh {

/// @brief Хранилище с Эйлеровыми / AMR-ячейками
using AmrStorage = Storage<zephyr::geom::AmrCell>;

}