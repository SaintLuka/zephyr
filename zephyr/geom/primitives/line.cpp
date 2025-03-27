#include <cstring>
#include <iostream>

#include <zephyr/geom/primitives/line.h>

namespace zephyr::geom {

inline double sqr(double x) {
    return x * x;
}

// Перпендикуляр единичной длины 'n' к вектору 'tau', лежащий в плоскости
// векторов 'a' и 'tau', сонаправленный с вектором 'a', то есть
// (n, n) = 1       // единичный
// (tau, n) = 0     // перпендикуляр
// (tau, n, a) = 0  // в одной плоскости
// (tau, a) > 0     // сонаправлены
inline Vector3d perpendicular(const Vector3d &tau, const Vector3d &a) {
    return tau.cross(a).cross(tau).normalized();
}

// ============================================================================
//                                    LINE
// ============================================================================

Line::Line(const Vector3d &v1, const Vector3d &v2)
        : verts{v1, v2} {}

Vector3d Line::operator()(double x) const {
    return get(verts[0], verts[1], x);
}

Vector3d Line::get(const Vector3d &v1, const Vector3d &v2, double x) {
    return 0.5 * (1.0 - x) * v1 + 0.5 * (1.0 + x) * v2;
}

Vector3d Line::normal(const Vector3d &v1, const Vector3d &v2, const Vector3d &c) {
    return perpendicular(v2 - v1, v1 - c);
}

double Line::Jacobian() const {
    return 0.5 * (verts[1] - verts[0]).norm();
}

Vector3d Line::center() const {
    return 0.5 * (verts[0] + verts[1]);
}

Vector3d Line::centroid() const {
    return 0.5 * (verts[0] + verts[1]);
}

// Точная формула, получена интегрированием
Vector3d Line::centroid(bool axial) const {
    Vector3d c = 0.5 * (verts[0] + verts[1]);
    if (!axial || c.y() == 0.0) { return c; }

    double dx = verts[1].x() - verts[0].x();
    double dy = verts[1].y() - verts[0].y();

    double coeff = dy / (12.0 * c.y());
    return Vector3d{c.x() + coeff * dx, c.y() + coeff * dy, c.z()};
}

Vector3d Line::normal(const Vector3d &c) const {
    return perpendicular(verts[1] - verts[0], verts[0] - c);
}

double Line::length() const {
    return (verts[1] - verts[0]).norm();
}

double Line::area_as() const {
    return 0.5 * (verts[0].y() + verts[1].y()) * (verts[1] - verts[0]).norm();
}

Vector3d Line::area_n(const Vector3d& c) const {
    return length() * normal(c);
}

// ============================================================================
//                                  SQ-LINE
// ============================================================================

SqLine::SqLine(const Vector3d &v1, const Vector3d &v2)
        : verts({v1, 0.5 * (v1 + v2), v2}) {}

SqLine::SqLine(const Vector3d &v1, const Vector3d &v2, const Vector3d &v3)
        : verts({v1, v2, v3}) {}

SqLine::SqLine(const Line &vs)
        : verts({vs[0], vs.center(), vs[1]}) {}

Vector3d SqLine::get(const Vector3d &v1, const Vector3d &vc, const Vector3d &v2, double x) {
    return vc + 0.5 * (v2 - v1) * x + 0.5 * (v2 - 2.0 * vc + v1) * x * x;
}

Vector3d SqLine::tangent(const Vector3d &v1, const Vector3d &vc, const Vector3d &v2, double x) {
    return 0.5 * (v2 - v1) + (v2 - 2.0 * vc + v1) * x;
}

Vector3d SqLine::projection(const Vector3d &a) const {
    const double eps = 1.0e-14;
    double L2 = (verts[2] - verts[0]).squaredNorm();

    // Перпендикуляр к плоскости, в которой лежит кривая
    // (ноль, если точки лежат на одной линии)
    Vector3d np = (verts[2] - verts[1]).cross(verts[0] - verts[1]);
    double np_norm = np.norm();

    if (np_norm < eps * L2) {
        // Считаем, что точки лежат на одной линии
        return a;
    } else {
        // Проекция точки 'a' на плоскость кривой
        np /= np_norm;
        return np.cross(a).cross(np);
    }
}

Vector3d SqLine::operator()(double x) const {
    return SqLine::get(verts[0], verts[1], verts[2], x);
}

Vector3d SqLine::get(double x) const {
    return SqLine::get(verts[0], verts[1], verts[2], x);
}

Vector3d SqLine::tangent(double x) const {
    return SqLine::tangent(verts[0], verts[1], verts[2], x);
}

Vector3d SqLine::normal(double x, const Vector3d &c) const {
    Vector3d a = projection(verts[1] - c);
    return perpendicular(tangent(x), a);
}

double SqLine::Jacobian(double x) const {
    return SqLine::tangent(verts[0], verts[1], verts[2], x).norm();
}

const Vector3d &SqLine::center() const {
    return verts[1];
}

Vector3d SqLine::normal(const Vector3d &c) const {
    // Это не заглушка, было установлено, что в расчете потока через
    // криволинейную грань следует использовать такую же нормаль,
    // как при расчете через обычную грань.
    // Единственное, здесь добавляется проекция на плоскость кривой.

    Vector3d a = projection(verts[1] - c);
    return perpendicular(verts[2] - verts[0], a);
}

double SqLine::length() const {
    // Это не заглушка, было установлено, что в расчете потока
    // по криволинейной грани следует использовать такую длину
    return (verts[2] - verts[0]).norm();
}

} // namespace zephyr::geom