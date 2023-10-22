#include <zephyr/mesh/mesh.h>
#include <zephyr/geom/grid.h>


namespace zephyr::mesh {


void Mesh::initialize(const Grid& grid) {
    m_locals.resize(grid.n_cells());

    for (int i = 0; i < grid.n_cells(); ++i) {
        m_locals[i].geom() = grid.amr_cell(i);
    }
}

} // namespace zephyr::mesh