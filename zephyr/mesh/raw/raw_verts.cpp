#include <zephyr/mesh/raw/raw_verts.h>
#include <zephyr/utils/mpi.h>

namespace zephyr::mesh {

RawVerts::RawVerts(bool unique_nodes)
    : has_nodes_(unique_nodes) {

}

void RawVerts::resize(index_t n_cells, index_t n_verts) {
    offsets.resize(n_cells + 1, offsets.back());

    coord.resize(n_verts);
    if (has_nodes_) {
        rank.resize(n_verts, utils::mpi::rank());
        index.resize(n_verts, -13);
        ghost.resize(n_verts, -1);
    }    
}

void RawVerts::resize_amr(index_t n_cells, int dim) {
    z_assert(dim == 2 || dim == 3, "RawVerts::resize_amr: bad dimension");
    int verts_per_cell = dim < 3 ? 9 : 27;
    resize(n_cells, verts_per_cell * n_cells);
}

void RawVerts::reserve(index_t n_cells, index_t n_verts) {
    offsets.reserve(n_cells + 1);

    coord.reserve(n_verts);
    if (has_nodes_) {
        rank.resize(n_verts);
        index.reserve(n_verts);
        ghost.reserve(n_verts);
    }    
}

void RawVerts::reserve_amr(index_t n_cells, int dim) {
    z_assert(dim == 2 || dim == 3, "RawVerts::resize_amr: bad dimension");
    int verts_per_cell = (dim < 3 ? 9 : 27);
    reserve(n_cells, verts_per_cell * n_cells);
}

void RawVerts::shrink_to_fit() {
    offsets.shrink_to_fit();

    coord.shrink_to_fit();
    if (has_nodes_) {
        rank.shrink_to_fit();
        index.shrink_to_fit();
        ghost.shrink_to_fit();
    }
}

memory_t RawVerts::memory_usage() const {
    memory_t mem;
    mem.add(offsets);
    mem.add(coord);
    mem.add(rank);
    mem.add(index);
    mem.add(ghost);
    return mem;
}

} // namespace zephyr::mesh
