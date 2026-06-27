#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/geom/surface/solid_body_2d.h>

namespace zephyr::geom {

// двумерная матрица поворота на угол phi
inline Matrix3d rotation_matrix(double phi) {
    double cos = std::cos(phi);
    double sin = std::sin(phi);
    Matrix3d R;
    R << cos, -sin, 0.0, sin, cos, 0.0, 0.0, 0.0, 1.0;
    return R;
}

SolidBody2D::SolidBody2D() { }

SolidBody2D::SolidBody2D(const Vector3d& center)
    : SolidBody(center) { }

void SolidBody2D::rotation_local(double angle) {
    SolidBody::rotation_local(rotation_matrix(angle));
}

void SolidBody2D::rotation_global(double angle) {
    SolidBody::rotation_global(rotation_matrix(angle));
}

void SolidBody2D::rotation_relative(double angle, const Vector3d& c) {
    SolidBody::rotation_relative(rotation_matrix(angle), c);
}

std::vector<Vector3d> SolidBody2D::outline(int n_points) const {
    auto vs = local_outline(n_points);
    for (Vector3d& v: vs) {
        v.applyOnTheLeft(m_rotation);
        v += m_center;
    }
    return vs;
}

// ============================================================================
//                                  STRIP
// ============================================================================

void BodyStrip::init_inside() {
    m_inside = [this](const Vector3d& v) -> bool {
        Vector3d p = this->in_local(v);
        return std::abs(p.x()) < 0.5 * m_width;
    };
}

BodyStrip::BodyStrip(double width) : m_width(width) {
    init_inside();
}

BodyStrip::BodyStrip(double width, const Vector3d& center)
    : SolidBody2D(center), m_width(width) {
    init_inside();
}

std::vector<Vector3d> BodyStrip::local_outline(int n_points) const {
    double a = 0.5 * m_width;
    double b = 4.0 * m_width;
    double c = 5.0 * m_width;
    return {
        Vector3d{-a, -b, 0.0},
        Vector3d{0., -c, 0.0},
        Vector3d{+a, -b, 0.0},
        Vector3d{+a, +b, 0.0},
        Vector3d{0., +c, 0.0},
        Vector3d{-a, +b, 0.0},
    };
}

// ============================================================================
//                                  DISK
// ============================================================================

void BodyDisk::init_inside() {
    m_inside = [this](const Vector3d &v) -> bool {
        return (v - m_center).norm() < m_radius;
    };
}

BodyDisk::BodyDisk(double radius) : m_radius(radius) {
    init_inside();
}

BodyDisk::BodyDisk(double radius, const Vector3d& center)
    : SolidBody2D(center), m_radius(radius) {
    init_inside();
}

std::vector<Vector3d> BodyDisk::local_outline(int n_points) const {
    std::vector<Vector3d> vs(n_points);
    for (int i = 0; i < n_points; ++i) {
        double phi = 2.0 * M_PI * (i - 1.0) / n_points;
        vs[i].x() = m_radius * std::cos(phi);
        vs[i].y() = m_radius * std::sin(phi);
        vs[i].z() = 0.0;
    }
    return vs;
}

double BodyDisk::volume_fraction(const mesh::EuCell& cell, double eps) const {
    auto poly = cell.polygon();
    return poly.disk_clip_area(m_center, m_radius) / cell.volume();
}

double BodyDisk::volume_inside(const mesh::EuCell& cell, double eps) const {
    auto poly = cell.polygon();
    return poly.disk_clip_area(m_center, m_radius);
}

// ============================================================================
//                                  SQUARE
// ============================================================================

void BodySquare::init_inside() {
    m_inside = [this](const Vector3d &v) -> bool {
        Vector3d p = this->in_local(v);
        return std::abs(p.x()) < 0.5 * m_length && std::abs(p.y()) < 0.5 * m_length;
    };
}

BodySquare::BodySquare(double length) : m_length(length) {
    init_inside();
}

BodySquare::BodySquare(double length, const Vector3d& center)
    : SolidBody2D(center), m_length(length) {
    init_inside();
}

std::vector<Vector3d> BodySquare::local_outline(int n_points) const {
    double a = 0.5 * m_length;
    return {
        Vector3d{-a, -a, 0.0},
        Vector3d{+a, -a, 0.0},
        Vector3d{+a, +a, 0.0},
        Vector3d{-a, +a, 0.0},
    };
}

} // namespace zephyr::geom