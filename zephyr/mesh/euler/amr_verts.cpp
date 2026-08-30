#include <zephyr/mesh/euler/amr_verts.h>
#include <zephyr/utils/mpi.h>

namespace zephyr::mesh {

void AmrVerts::resize(index_t n_verts) {
    coord.resize(n_verts);
    if (unique_) {
        rank.resize(n_verts);
        index.resize(n_verts);
        ghost.resize(n_verts);
    }    
}

void AmrVerts::reserve(index_t n_verts) {
    coord.reserve(n_verts);
    if (unique_) {
        rank.resize(n_verts);
        index.reserve(n_verts);
        ghost.reserve(n_verts);
    }    
}

void AmrVerts::shrink_to_fit() {
    coord.shrink_to_fit();
    if (!unique_) {
        rank.clear();
        index.clear();
        ghost.clear();
    }
    rank.shrink_to_fit();
    index.shrink_to_fit();
    ghost.shrink_to_fit();
}

void AmrVerts::clear_unique() {
    unique_ = false;
    rank.clear();
    index.clear();
    ghost.clear();
}

void AmrVerts::init_unique(index_t idx, index_t gst) {
    unique_ = true;
    rank.clear();
    rank.resize(coord.size(), utils::mpi::rank());
    index.clear();
    index.resize(coord.size(), idx);
    ghost.clear();
    ghost.resize(coord.size(), gst);
}

memory_t AmrVerts::memory_usage() const {
    memory_t mem;
    mem.add(offsets);
    mem.add(coord);
    mem.add(rank);
    mem.add(index);
    mem.add(ghost);
    return mem;
}

} // namespace zephyr::mesh
