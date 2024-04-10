#include <zephyr/geom/box.h>


namespace zephyr { namespace geom {


Box::Box(const Vector3d &_vmin, const Vector3d &_vmax)
        : vmin(_vmin), vmax(_vmax) {}


Vector3d Box::center() const {
    return 0.5 * (vmin + vmax);
}

Vector3d Box::size() const {
    return vmax - vmin;
}

double Box::diameter() const {
    return (vmax - vmin).norm();
}

bool Box::inside(const Vector3d& p) const {
    return vmin.x() <= p.x() && p.x() <= vmax.x() &&
           vmin.y() <= p.y() && p.y() <= vmax.y() &&
           vmin.z() <= p.z() && p.z() <= vmax.z();
}

void Box::extend(double margin) {
    extend(margin, margin, margin);
}

void Box::extend(double margin_x, double margin_y, double margin_z) {
    Vector3d L = size();
    vmin.x() -= margin_x * L.x();
    vmin.y() -= margin_y * L.y();
    vmin.z() -= margin_z * L.z();

    vmax.x() += margin_x * L.x();
    vmax.y() += margin_y * L.y();
    vmax.z() += margin_z * L.z();
}

QuasiRandom2D Box::quasiRandom2D() const {
    return QuasiRandom2D(vmin, size());
}

QuasiRandom2D::QuasiRandom2D(const Vector3d &_vmin, const Vector3d &_size)
    : vmin(_vmin), size(_size) {

    const double phi2 = 1.32471795724474602596;
    step = {size.x() / phi2, size.y() / (phi2 * phi2), 0.0};

    vmin.z() = 0.0;
    size.z() = 0.0;

    shift = Vector3d::Zero();
}

Vector3d QuasiRandom2D::get() {
    Vector3d res = vmin + shift;

    shift += step;
    shift.x() = std::fmod(shift.x(), size.x());
    shift.y() = std::fmod(shift.y(), size.y());

    return res;
}

} // zephyr
} // geom