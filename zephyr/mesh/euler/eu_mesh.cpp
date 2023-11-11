
#include <zephyr/geom/grid.h>
#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/geom/primitives/amr_faces.h>
#include <zephyr/mesh/euler/eu_cell.h>
#include <zephyr/mesh/euler/eu_mesh.h>

namespace zephyr::mesh {


void EuMesh::initialize(const Grid& grid) {
    m_locals.resize(grid.n_cells());

    for (int i = 0; i < grid.n_cells(); ++i) {
        m_locals[i] = grid.amr_cell(i);
    }
}

} // namespace zephyr::mesh