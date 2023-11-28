
#include <zephyr/geom/grid.h>
#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/geom/primitives/bfaces.h>
#include <zephyr/geom/box.h>
#include <zephyr/mesh/euler/eu_cell.h>
#include <zephyr/mesh/euler/eu_mesh.h>

namespace zephyr::mesh {


void EuMesh::initialize(const Grid& grid) {
    m_locals.resize(grid.n_cells());

    for (int i = 0; i < grid.n_cells(); ++i) {
        m_locals[i] = grid.amr_cell(i);
    }
}

geom::Box EuMesh::bbox() {
    double inf = std::numeric_limits<double>::infinity();
    Vector3d vmin = {+inf, +inf, +inf};
    Vector3d vmax = {-inf, -inf, -inf};

    for (auto& cell: *this) {
        for (auto& face: cell.faces()) {
            for (int i = 0; i < face.size(); ++i) {
                auto& v = face.vs(i);
                vmin = vmin.cwiseMin(v);
                vmax = vmax.cwiseMax(v);
            }
        }
    }

    return geom::Box(vmin, vmax);
}

} // namespace zephyr::mesh