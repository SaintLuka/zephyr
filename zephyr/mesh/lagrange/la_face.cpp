#include <cstring>

#include <zephyr/geom/primitives/mov_cell.h>
#include <zephyr/geom/primitives/bfaces.h>
#include <zephyr/mesh/lagrange/la_face.h>
#include <zephyr/mesh/lagrange/la_cell.h>

namespace zephyr::mesh {

LaFace::LaFace(const LaCell &cell, BFace* self, BFace* end, Direction dir)
        : m_cell(cell), m_face(self), m_end(end), m_dir(dir) {

    // Ищем первую определенную грань
    while (m_face < end && m_face->to_skip(m_dir)) {
        ++m_face;
    }
}

const Vector3d &LaFace::normal() const {
    return m_face->normal;
}

LaCell LaFace::neib() const {
    return m_face->is_boundary() ? m_cell : m_cell.neib(*m_face);
}

const Byte* LaFace::neib_data() const {
    return m_face->is_boundary() ? m_cell.data() : m_cell.neib(*m_face).data();
}

Boundary LaFace::flag() const {
    return m_face->boundary;
}

bool LaFace::is_boundary() const {
    return m_face->is_boundary();
}

void LaFace::set_boundary(Boundary flag) {
    m_face->boundary = flag;
}

double LaFace::area() const {
    return m_face->area;
}

Vector3d LaFace::center() const {
    return m_face->center;
}
/*
const Vector3d& LaFace::vs(int idx) const {
    return m_cell.vs(m_face->vertices[idx]);
}
*/
LaFaces::LaFaces(const LaCell& cell, Direction dir)
    : m_cell(cell), m_dir(dir) {

}

LaFace LaFaces::begin() const {
    auto &cell = const_cast<CellStorage::Iterator &>(m_cell.m_it);
    BFace *beg = cell->faces.begin();
    BFace *end = cell->faces.end();
    return {m_cell, beg, end, m_dir};
}

LaFace LaFaces::end() const {
    auto &cell = const_cast<CellStorage::Iterator &>(m_cell.m_it);
    BFace *end = cell->faces.end();
    return {m_cell, end, end, m_dir};
}

} // namespace zephyr::mesh
