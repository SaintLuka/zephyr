#include <zephyr/mesh/euler/amr_verts.h>

namespace zephyr::mesh {

void AmrVerts::resize(index_t n_verts) {
    coords.resize(n_verts);
    if (m_unique) {
        index.resize(n_verts);
        ghost.resize(n_verts);
    }    
}

void AmrVerts::reserve(index_t n_verts) {
    coords.reserve(n_verts);
    if (m_unique) {
        index.reserve(n_verts);
        ghost.reserve(n_verts);
    }    
}

void AmrVerts::shrink_to_fit() {
    coords.shrink_to_fit();
    if (!m_unique) {
        index.clear();
        ghost.clear();
    }
    index.shrink_to_fit();
    ghost.shrink_to_fit();
}

void AmrVerts::clear_unique() {
    m_unique = false;
    index.clear();
    ghost.clear();
}

void AmrVerts::init_unique() {
    m_unique = true;
    index.clear();
    index.resize(coords.size(), -1);
    ghost.clear();
    ghost.resize(coords.size(), -1);
}

memory_t AmrVerts::memory_usage() const {
    memory_t mem;
    mem.add(coords);
    mem.add(index);
    mem.add(ghost);
    return mem;
}

} // namespace zephyr::mesh
