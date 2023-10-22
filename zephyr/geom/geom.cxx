#include <zephyr/geom/geom.h>

namespace zephyr::geom {

// ============================================================================
//                                  POLYGON
// ============================================================================


double area(const VerticesList &vs) {
    int n_points = vs.size();
    Vector3d cell_c = center(vs);

    double S = 0.0;
    for (int i = 0; i < n_points; ++i) {
        const Vector3d &v1 = vs[i];
        const Vector3d &v2 = vs[(i + 1) % n_points];

        Vector3d face_c = 0.5 * (v1 + v2);
        Vector3d normal = {v2.y() - v1.y(), v1.x() - v2.x(), 0.0};

        S += std::abs((face_c - cell_c).dot(normal));
    }
    S /= 2.0;

    return S;
}


Vector3d center(const VerticesList &vs) {
    Vector3d C = {0.0, 0.0, 0.0};
    for (auto &v: vs) {
        C += v;
    }
    return C / vs.size();
}

Vector3d centroid(const VerticesList &vs, double area) {
    if (area == 0.0) {
        area = geom::area(vs);
    }
    int n_points = vs.size();

    Vector3d C = {0.0, 0.0, 0.0};
    for (int i = 0; i < n_points; ++i) {
        auto &v1 = vs[i];
        auto &v2 = vs[(i + 1) % n_points];

        C.x() += (v2.y() - v1.y()) * (v2.x() * v2.x() + v2.x() * v1.x() + v1.x() * v1.x());
        C.y() -= (v2.x() - v1.x()) * (v2.y() * v2.y() + v2.y() * v1.y() + v1.y() * v1.y());
    }
    C /= (6.0 * area);

    return C;
}

} // namespace zephyr::geom