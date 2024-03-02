#pragma once

#include <functional>
#include <zephyr/geom/vector.h>
#include <zephyr/mesh/mesh.h>

namespace zephyr::math {

template<class T>
std::array<T, 3> compute_grad(Cell &cell, const std::function<T(Cell &)> &to_state) {
    std::array<T, 3> z_xyz;

    T zc = to_state(cell);

    for (auto &face: cell.faces()) {
        auto neib = face.neib();
        T zn = to_state(neib);

        Eigen::Vector3d S = face.normal() * face.area();

        double d1 = abs(1 / (face.center() - cell.center()).dot(face.normal()));
        double d2 = abs(1 / (face.center() - neib.center()).dot(face.normal()));
        double a1 = d1 / (d1 + d2), a2 = d2 / (d1 + d2);
        T state = a1 * zc.vec() + a2 * zn.vec();

        z_xyz[0].vec() += state.vec() * S.x();
        z_xyz[1].vec() += state.vec() * S.y();
        z_xyz[2].vec() += state.vec() * S.z();
    }

    for (auto &z: z_xyz)
        z.vec() /= cell.volume();

    return z_xyz;
}

}