#include <zephyr/geom/primitives/bvertices.h>

namespace zephyr::geom {

BVertices::BVertices(const Polygon &poly) {
    int mid_count = std::min(poly.size(), max_count);

    for (int i = 0; i < mid_count; ++i) {
        verts[i] = poly[i];
    }
    for (int i = mid_count; i < max_count; ++i) {
        verts[i] = {NAN, NAN, NAN};
    }
};

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

} // namespace zephyr::geom