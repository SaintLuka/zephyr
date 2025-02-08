#include <zephyr/geom/box.h>


namespace zephyr::geom {

Box::Box()
    : vmin{NAN, NAN, NAN},
      vmax{NAN, NAN, NAN} {

}

Box::Box(const Vector3d &_vmin, const Vector3d &_vmax)
    : vmin(_vmin), vmax(_vmax) {

}

Vector3d Box::center() const {
    return 0.5 * (vmin + vmax);
}

Vector3d Box::size() const {
    return vmax - vmin;
}

bool Box::is_2D() const {
    return vmax.z() == vmin.z();
}

bool Box::is_3D() const {
    return vmax.z() > vmin.z();
}

double Box::diameter() const {
    return (vmax - vmin).norm();
}

double Box::area() const {
    return (vmax[0] - vmin[0]) * (vmax[1] - vmin[1]);
}

double Box::volume() const {
    return (vmax[0] - vmin[0]) * (vmax[1] - vmin[1]) * (vmax[2] - vmin[2]);
}

bool Box::inside(const Vector3d& p) const {
    return vmin.x() <= p.x() && p.x() <= vmax.x() &&
           vmin.y() <= p.y() && p.y() <= vmax.y() &&
           vmin.z() <= p.z() && p.z() <= vmax.z();
}

Vector3d Box::shove_in(const Vector3d& p) const {
    return {
            std::max(vmin.x(), std::min(p.x(), vmax.x())),
            std::max(vmin.y(), std::min(p.y(), vmax.y())),
            std::max(vmin.z(), std::min(p.z(), vmax.y()))
    };
}

std::vector<double> Box::outline_x() const {
    return {vmin.x(), vmax.x(), vmax.x(), vmin.x(), vmin.x()};
}

std::vector<double> Box::outline_y() const {
    return {vmin.y(), vmin.y(), vmax.y(), vmax.y(), vmin.y()};
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

Random2D Box::random2D(int seed) const {
    return Random2D(vmin, vmax, seed);
}

QuasiRandom2D Box::quasiRandom2D() const {
    return QuasiRandom2D(vmin, size());
}

std::ostream& operator<<(std::ostream& os, const Box& box) {
    os << "[" << box.vmax.x() << ", " << box.vmax.x() << "] × [" <<
       box.vmax.x() << ", " << box.vmax.x() << "]";
    if (box.vmax.z() - box.vmin.z() > 0.0) {
        os << " × [" << box.vmax.z() << ", " << box.vmax.z() << "]";
    }
    return os;
}


Random2D::Random2D(const Vector3d& vmin, const Vector3d& vmax, int seed) {
    gen = std::mt19937_64(seed);

    distr_x = std::uniform_real_distribution<double>(vmin.x(), vmax.x());
    distr_y = std::uniform_real_distribution<double>(vmin.y(), vmax.y());
}

Vector3d Random2D::get() {
    return {distr_x(gen), distr_y(gen), 0.0};
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

} // namespace zephyr::geom