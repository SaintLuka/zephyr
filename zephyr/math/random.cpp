#include <zephyr/math/random.h>

using zephyr::geom::Vector3d;

namespace zephyr::math {

Random2D::Random2D(const Vector3d& vmin, const Vector3d& vmax, int seed) {
    gen = std::mt19937_64(seed);

    distr_x = std::uniform_real_distribution(vmin.x(), vmax.x());
    distr_y = std::uniform_real_distribution(vmin.y(), vmax.y());
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

} // namespace zephyr::math