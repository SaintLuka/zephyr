#include <zephyr/mesh/primitives/bvertices.h>

using namespace zephyr::geom;

namespace zephyr::mesh {

BVertices::BVertices(const Polygon &poly) {
    int mid_count = std::min(poly.size(), max_count);

    for (int i = 0; i < mid_count; ++i) {
        verts[i] = poly[i];
    }
    for (int i = mid_count; i < max_count; ++i) {
        verts[i] = {NAN, NAN, NAN};
    }
}

BVertices::BVertices(const Polyhedron &poly) {
    int mid_count = std::min(poly.n_verts(), max_count);

    for (int i = 0; i < mid_count; ++i) {
        verts[i] = poly.vertex(i);
    }
    for (int i = mid_count; i < max_count; ++i) {
        verts[i] = {NAN, NAN, NAN};
    }
}

int BVertices::find(const Vector3d &v1, double eps) const {
    double eps2 = eps * eps;

    for (int iv = 0; iv < BVertices::max_count; ++iv) {
        const Vector3d& v2 = verts[iv];

        // Нашли интересующую нас вершину
        if (!v2.hasNaN() && (v2 - v1).squaredNorm() < eps2) {
            return iv;
        }
    }

    return -1;
}

} // namespace zephyr::mesh