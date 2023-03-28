#include <zephyr/mesh/generator/vertex.h>
#include <zephyr/mesh/generator/curve/circle.h>

namespace zephyr { namespace mesh { namespace generator {

Circle::Circle(double _R, const Vector3d &_v0) {
    R = _R;
    v0 = _v0;
}

Curve::Ptr Circle::create(double _R, const Vector3d &_v0) {
    return std::make_shared<Circle>(_R, _v0);
}

Curve::Ptr Circle::create(BaseVertex::Ref v1, BaseVertex::Ref v2, BaseVertex::Ref v3) {
    Vector3d a = v2->v() - v1->v();
    Vector3d b = v3->v() - v1->v();

    double L = std::max(a.norm(), b.norm());

    double cross = std::abs(a.cross(b).z());

    if (std::sqrt(cross) < 1.0e-10 * L) {
        throw std::runtime_error("Can't create circle (three points belong to the same line)");
    }

    Vector3d B = v2->v() - v1->v();
    Vector3d C = v3->v() - v1->v();

    double D = 2.0 * (B.cross(C).z());

    Vector3d v0 = {
        v1->v().x() + (C.y() * (B.x() * B.x() + B.y() * B.y()) - B.y() * (C.x() * C.x() + C.y() * C.y())) / D,
        v1->v().y() + (B.x() * (C.x() * C.x() + C.y() * C.y()) - C.x() * (B.x() * B.x() + B.y() * B.y())) / D,
        0.0
    };

    double R = ((v0 - v1->v()).norm() + (v0 - v2->v()).norm() + (v0 - v3->v()).norm()) / 3.0;

    double err =
        std::fabs((v0 - v1->v()).norm() / R - 1.0) +
        std::fabs((v0 - v2->v()).norm() / R - 1.0) +
        std::fabs((v0 - v3->v()).norm() / R - 1.0);

    if (err > 1.0e-10) {
        throw std::runtime_error("Can't create circle");
    }

    return std::make_shared<Circle>(R, v0);
}

double Circle::get_X(double t) const {
    double phi = M_PI * t;
    return v0.x() + R * std::cos(phi);
}

double Circle::get_Y(double t) const {
    double phi = M_PI * t;
    return v0.y() + R * std::sin(phi);
}

Vector3d Circle::get_V(double t) const {
    double phi = M_PI * t;
    return {
        v0.x() + R * std::cos(phi),
        v0.y() + R * std::sin(phi),
        0.0
    };
}

Vector3d Circle::get_T(double t) const {
    Vector3d v = get_V(t);
    Vector3d T = {v.y() - v0.y(), -v.x() + v0.x(), 0.0};
    return T.normalized();
}

Vector3d Circle::get_N(double t) const {
    Vector3d v = get_V(t);
    Vector3d N = {v.x() - v0.x(), v.y() - v0.y(), 0.0};
    return N.normalized();
}

Vector3d Circle::projection(const Vector3d &v) const {
    Vector3d normal = v - v0;
    normal.normalize();
    normal *= R;
    return v0 + normal;
}

} // namespace generator
} // namespace mesh
} // namespace zephyr
