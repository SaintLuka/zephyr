#include <zephyr/mesh/euler/soa_prim.h>
#include <zephyr/geom/geom.h>

namespace zephyr::mesh {

using namespace geom;



geom::Polygon QCell::polygon() const {
    return m_cells->polygon(m_index);
}

} //