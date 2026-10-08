#include <zephyr/mesh/cell.h>

#include <zephyr/geom/surface/solid_body.h>

namespace zephyr::geom {

inline bool basic_inside(const Vector3d& ) {
    return false;
}

SolidBody::SolidBody()  :
    m_center(Vector3d::Zero()),
    m_rotation(Matrix3d::Identity()),
    m_inside(basic_inside) {

}

SolidBody::SolidBody(const Vector3d& center) :
    m_center(center),
    m_rotation(Matrix3d::Identity()),
    m_inside(basic_inside) {

}

Vector3d SolidBody::center() const {
    return m_center;
}

void SolidBody::move(const Vector3d& shift) {
    m_center += shift;
}

void SolidBody::rotation_local(const Matrix3d& R) {
    m_rotation.applyOnTheLeft(R);
}

void SolidBody::rotation_global(const Matrix3d& R) {
    m_rotation.applyOnTheLeft(R);
    m_center.applyOnTheLeft(R);
}

void SolidBody::rotation_relative(const Matrix3d& R, const Vector3d& c) {
    m_rotation.applyOnTheLeft(R);
    m_center -= c;
    m_center.applyOnTheLeft(R);
    m_center += c;
}

bool SolidBody::inside(const Vector3d& v) const {
    return m_inside(v);
}

double SolidBody::volume_fraction(const mesh::Cell& cell, double eps) const {
    double vol_frac = cell.approx_vol_fraction(m_inside);
    if (0.0 < vol_frac && vol_frac < 1.0) {
        int n_points = std::max(4, std::min(static_cast<int>(1.0 / eps), 100'000'000));
        vol_frac = cell.volume_fraction(m_inside, n_points);
    }
    return vol_frac;
}

double SolidBody::volume_inside(const mesh::Cell& cell, double eps) const {
    return volume_fraction(cell, eps) * cell.volume();
}

Vector3d SolidBody::in_local(const Vector3d& v) {
    return m_rotation.transpose() * (v - m_center);
}

} // namespace zephyr::geom