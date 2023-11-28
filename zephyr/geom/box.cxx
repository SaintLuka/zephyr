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

} // zephyr
} // geom