#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/geom/geom.h>

namespace zephyr::mesh {

using namespace geom;



geom::Polygon EuCell::polygon() const {
    return m_cells->polygon(m_index);
}

} //