#include <zephyr/mesh/primitives/Side3D.h>
#include <zephyr/mesh/primitives/mov_cell.h>
#include <zephyr/mesh/lagrange/la_cell.h>

namespace zephyr::mesh {

CellStorage::Iterator safe_iterator(CellStorage &locals, int idx) {
    if (idx < locals.size()) {
        return locals.iterator(idx);
    } else {
        return locals.end();
    }
}

CellStorage::Iterator safe_iterator(CellStorage &locals, const Adjacent &adj) {
    if (adj.index < locals.size()) {
        return locals.iterator(adj.index);
    } else {
        return locals.end();
    }
}

LaCell::LaCell(CellStorage &_locals, int idx)
    : m_it(safe_iterator(_locals, idx)), m_locals(_locals) {

}

LaCell::LaCell(CellStorage &_locals, const Adjacent &adj)
    : m_it(safe_iterator(_locals, adj)), m_locals(_locals) {

}

LaCell LaCell::locals(int idx) {
    return { m_locals, idx };
}

LaCell LaCell::neib(const BFace &face) const {
    return { m_locals, face.adjacent };
}

LaCell LaCell::neib(const LaFace &face) const {
    return { m_locals, face.adjacent() };
}

LaCell LaCell::neib(int face_idx) const {
    return { m_locals, m_it->faces[face_idx].adjacent };
}

LaFaces LaCell::faces(Direction dir) const {
    return { *this, dir };
}

} // namespace zephyr::mesh
