/// @file Простое, но очень полезное определеине.

#include <zephyr/mesh/primitives/amr_cell.h>
#include <zephyr/mesh/storage.h>

namespace zephyr::mesh {

/// @brief Хранилище с Эйлеровыми / AMR-ячейками
using AmrStorage = Storage<AmrCell>;

}