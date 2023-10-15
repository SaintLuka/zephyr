#include <zephyr/geom/generator/curve/curve.h>

namespace zephyr::geom::generator {

Curve::Curve()
    : m_flag(Boundary::UNDEFINED) {
}

Vector3d Curve::get_V(double t) const {
    return Vector3d(get_X(t), get_Y(t), 0.0);
}

Vector3d Curve::get_T(double t) const {
    const double dt = 1.0e-8;
    Vector3d T = {
        get_X(t + dt) - get_X(t - dt),
        get_Y(t + dt) - get_Y(t - dt),
        0.0
    };
    T.normalize();
    return T;
}

Vector3d Curve::get_N(double t) const {
    const double dt = 1.0e-8;
    Vector3d N = {
        get_Y(t + dt) - get_Y(t - dt),
        get_X(t - dt) - get_X(t + dt),
        0.0
    };
    N.normalize();
    return N;
}

void Curve::set_boundary_flag(Boundary flag) {
    m_flag = flag;
}

Boundary Curve::boundary_flag() const {
    return m_flag;
}

} // namespace zephyr::geom::generator
