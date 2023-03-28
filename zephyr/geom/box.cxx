#include <zephyr/geom/box.h>


namespace zephyr { namespace geom {


Box::Box(const Vector3d &_vmin, const Vector3d &_vmax)
        : vmin(_vmin), vmax(_vmax) {}


} // zephyr
} // geom