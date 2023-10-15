#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/line.h>

namespace zephyr::geom::generator {

Line::Line(const Vector3d &_v1, const Vector3d &_v2)
    : v1(_v1), v2(_v2) {
}

Curve::Ptr Line::create(const Vector3d &_v1, const Vector3d &_v2) {
    return std::make_shared<Line>(_v1, _v2);
}

Curve::Ptr Line::create(BaseVertex::Ref v1, BaseVertex::Ref v2) {
    return std::make_shared<Line>(v1->v(), v2->v());
}

double Line::get_X(double t) const {
    double xi = std::tan(M_PI_2 * t);
    return 0.5 * (v1.x() + v2.x()) + 0.5 * (v2.x() - v1.x()) * xi;
}

double Line::get_Y(double t) const {
    double xi = std::tan(M_PI_2 * t);
    return 0.5 * (v1.y() + v2.y()) + 0.5 * (v2.y() - v1.y()) * xi;
}

Vector3d Line::get_V(double t) const {
    double xi = std::tan(M_PI_2 * t);
    return {
        0.5 * (v1.x() + v2.x()) + 0.5 * (v2.x() - v1.x()) * xi,
        0.5 * (v1.y() + v2.y()) + 0.5 * (v2.y() - v1.y()) * xi,
        0.0
    };
}

Vector3d Line::get_T(double t) const {
    Vector3d T = {v2.x() - v1.x(), v2.y() - v1.y(), 0.0};
    return T.normalized();
}

Vector3d Line::get_N(double t) const {
    Vector3d N = {v2.y() - v1.y(), v1.x() - v2.x(), 0.0};
    return N.normalized();
}

Vector3d Line::projection(const Vector3d &v) const {
    Vector3d n = (v2 - v1).normalized();
    return v1 + n * (v - v1).dot(n);
}

} // namespace zephyr::geom::generator
