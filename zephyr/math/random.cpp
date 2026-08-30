#include <zephyr/math/random.h>

using zephyr::geom::Vector3d;

namespace zephyr::math {

Random::Ptr Random::Uniform(const geom::Box& box, int seed) {
    if (box.is_2D()) {
        return std::make_shared<Random2D>(box.vmin, box.vmax, seed);
    }
    else {
        return std::make_shared<Random3D>(box.vmin, box.vmax, seed);
    }
}

Random::Ptr Random::Quasi(const geom::Box& box) {
    if (box.is_2D()) {
        return std::make_shared<QuasiRandom2D>(box.vmin, box.sizes());
    }
    else {
        throw std::runtime_error("Not implemented");
    }
}

Random2D::Random2D(const geom::Box &box, int seed)
    : Random2D(box.vmin, box.vmax, seed) { }

Random2D::Random2D(const Vector3d& vmin, const Vector3d& vmax, int seed) {
    gen = std::mt19937_64(seed);

    distr_x = std::uniform_real_distribution(vmin.x(), vmax.x());
    distr_y = std::uniform_real_distribution(vmin.y(), vmax.y());
}

Vector3d Random2D::next() {
    return {distr_x(gen), distr_y(gen), 0.0};
}

Random3D::Random3D(const geom::Box &box, int seed)
    : Random3D(box.vmin, box.vmax, seed) { }

Random3D::Random3D(const Vector3d& vmin, const Vector3d& vmax, int seed) {
    gen = std::mt19937_64(seed);

    distr_x = std::uniform_real_distribution(vmin.x(), vmax.x());
    distr_y = std::uniform_real_distribution(vmin.y(), vmax.y());
    distr_z = std::uniform_real_distribution(vmin.z(), vmax.z());
}

Vector3d Random3D::next() {
    return {distr_x(gen), distr_y(gen), distr_z(gen)};
}

QuasiRandom2D::QuasiRandom2D(const geom::Box& box)
    : QuasiRandom2D(box.vmin, box.sizes()) { }

QuasiRandom2D::QuasiRandom2D(const Vector3d &_vmin, const Vector3d &_size)
    : vmin(_vmin), size(_size) {

    const double phi2 = 1.32471795724474602596;
    step = {size.x() / phi2, size.y() / (phi2 * phi2), 0.0};

    vmin.z() = 0.0;
    size.z() = 0.0;

    shift = Vector3d::Zero();
}

Vector3d QuasiRandom2D::next() {
    Vector3d res = vmin + shift;

    shift += step;
    shift.x() = std::fmod(shift.x(), size.x());
    shift.y() = std::fmod(shift.y(), size.y());

    return res;
}

} // namespace zephyr::math