#include <cstring>

#include <zephyr/mesh/primitives/amr_cell.h>
#include <zephyr/mesh/primitives/bfaces.h>
#include <zephyr/mesh/euler/eu_face.h>
#include <zephyr/mesh/euler/eu_cell.h>

using namespace zephyr::geom;

namespace zephyr::mesh {

EuFace::EuFace(const EuCell &cell, Side side)
    : m_cell(cell),
      m_face(nullptr),
      m_end (nullptr),
      m_dir(Direction::ANY) {

    m_face = m_end = &const_cast<AmrCell &>(cell.geom()).faces[side];
}

EuFace::EuFace(const EuCell &cell, BFace* self, BFace* end, Direction dir)
        : m_cell(cell), m_face(self), m_end(end), m_dir(dir) {

    // Ищем первую определенную грань
    while (m_face < end && m_face->to_skip(m_dir)) {
        ++m_face;
    }
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

const Vector3d& EuFace::normal() const {
    return m_face->normal;
}

double EuFace::area() const {
    return m_face->area;
}

double EuFace::area(bool axial) const {
    return m_face->get_area(axial);
}

Vector3d EuFace::area_n() const {
    return m_face->area_n();
}

const Vector3d& EuFace::vs(int idx) const {
    return m_cell.vs(m_face->vertices[idx]);
}

Vector3d EuFace::symm_point(const Vector3d& p) const {
    return p + 2.0 * (m_face->center - p).dot(m_face->normal) * m_face->normal;
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
