#pragma once

#include <zephyr/geom/primitives/basic_face.h>

namespace zephyr::geom {

/// @brief AmrFace содержит не более 4 вершин
using AmrFace = BasicFace<4>;

} // namespace zephyr::geom