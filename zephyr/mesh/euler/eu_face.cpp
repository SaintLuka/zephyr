#include <cstring>

#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/geom/primitives/bfaces.h>
#include <zephyr/mesh/euler/eu_face.h>
#include <zephyr/mesh/euler/eu_cell.h>

namespace zephyr::mesh {

EuFace::EuFace(const EuCell &cell, geom::Side side)
    : m_cell(cell),
      m_face(nullptr),
      m_end (nullptr),
      m_dir(Direction::ANY) {

    m_face = m_end = &const_cast<geom::AmrCell &>(cell.geom()).faces[side];
}

EuFace::EuFace(const EuCell &cell, BFace* self, BFace* end, Direction dir)
        : m_cell(cell), m_face(self), m_end(end), m_dir(dir) {

    // Ищем первую определенную грань
    while (m_face < end && m_face->to_skip(m_dir)) {
        ++m_face;
    }
}

const Vector3d &EuFace::normal() const {
    return m_face->normal;
}

EuCell EuFace::neib() const {
    return m_face->is_boundary() ? m_cell : m_cell.neib(*m_face);
}

const Byte* EuFace::neib_data() const {
    return m_cell.neib_data(*m_face);
}

Boundary EuFace::flag() const {
    return m_face->boundary;
}

bool EuFace::is_boundary() const {
    return m_face->is_boundary();
}

void EuFace::set_boundary(Boundary flag) {
    m_face->boundary = flag;
}

double EuFace::area() const {
    return m_face->area;
}

const Vector3d& EuFace::vs(int idx) const {
    return m_cell.vs(m_face->vertices[idx]);
}

EuFaces::EuFaces(const EuCell& cell, Direction dir)
    : m_cell(cell), m_dir(dir) {

}

EuFace EuFaces::begin() const {
    auto &cell = const_cast<AmrStorage::Iterator &>(m_cell.m_it);
    BFace *beg = cell->faces.begin();
    BFace *end = cell->faces.end();
    return {m_cell, beg, end, m_dir};
}

EuFace EuFaces::end() const {
    auto &cell = const_cast<AmrStorage::Iterator &>(m_cell.m_it);
    BFace *end = cell->faces.end();
    return {m_cell, end, end, m_dir};
}

} // namespace zephyr::mesh
