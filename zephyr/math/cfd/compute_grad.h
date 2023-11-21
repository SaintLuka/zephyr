#pragma once

#include <functional>
#include <zephyr/geom/vector.h>
#include <zephyr/mesh/mesh.h>

namespace zephyr::math {

template<class T>
std::array<T, 3> compute_grad(ICell &cell, const std::function<T(zephyr::mesh::ICell &)> &to_state) {
    std::array<T, 3> z_xyz;

    T zc = to_state(cell);

    for (auto &face: cell.faces()) {
        auto neib = face.neib();
        T zn = to_state(neib);

        Eigen::Vector3d S = 0.5 * face.normal() * face.area();

        T state = zc.vec() + zn.vec();
        z_xyz[0].vec() += state.vec() * S.x();
        z_xyz[1].vec() += state.vec() * S.y();
        z_xyz[2].vec() += state.vec() * S.z();
    }

    for (auto &z: z_xyz)
        z.vec() /= cell.volume();

    return z_xyz;
}

}