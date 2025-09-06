#include <zephyr/phys/tests/test_3D.h>

namespace zephyr::phys {

// ============================================================================
//                         Sedov Blast 3D Test
// ============================================================================

SedovBlast3D::SedovBlast3D(Sedov3D::params p)
    : exact(p) {
    m_materials += IdealGas::create(p.gamma);
    
    init_time = exact.time_by_radius(0.4);
    finish = exact.time_by_radius(0.9);
}

double SedovBlast3D::density(const Vector3d &r) const {
    return density_t(r, init_time);    
}

Vector3d SedovBlast3D::velocity(const Vector3d &r) const {
    return velocity_t(r, init_time);    
}

double SedovBlast3D::pressure(const Vector3d &r) const {
    return pressure_t(r, init_time);
}

double SedovBlast3D::density_t(const Vector3d &r, double t) const {
    return exact.density(r.norm(), t);
}

Vector3d SedovBlast3D::velocity_t(const Vector3d &r, double t) const {
    return exact.velocity(r.norm(), t) * r.normalized();
}

double SedovBlast3D::pressure_t(const Vector3d &r, double t) const {
    return exact.pressure(r.norm(), t);    
}

} // namespace zephyr::phys