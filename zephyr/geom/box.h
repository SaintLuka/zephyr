#pragma once

#include <zephyr/geom/vector.h>

namespace zephyr::geom {

// TODO: Написать нормально
class Box {
public:

    Vector3d vmin, vmax;

    Box(const Vector3d& _vmin, const Vector3d& _vmax);
};

} // namespace zephyr::geom

