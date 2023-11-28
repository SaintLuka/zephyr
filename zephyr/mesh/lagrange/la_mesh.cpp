#include <zephyr/geom/grid.h>
#include <zephyr/geom/box.h>
#include <zephyr/mesh/lagrange/la_mesh.h>

namespace zephyr::mesh {

void LaMesh::initialize(const Grid& grid) {
    m_nodes.resize(grid.n_nodes());
    for (int i = 0; i < grid.n_nodes(); ++i) {
        m_nodes[i] = *grid.node(i);
    }

    m_cells.resize(grid.n_cells());
    for (int i = 0; i < grid.n_cells(); ++i) {
        m_cells[i] = grid.mov_cell(i);
    }
}

geom::Box LaMesh::bbox() {
    double inf = std::numeric_limits<double>::infinity();
    Vector3d vmin = {+inf, +inf, +inf};
    Vector3d vmax = {-inf, -inf, -inf};

    for (auto& node: m_nodes) {
        vmin = vmin.cwiseMin(node.coords);
        vmax = vmax.cwiseMax(node.coords);
    }

    return geom::Box(vmin, vmax);
}

} // namespace zephyr::mesh