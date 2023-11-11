#include <cstring>

#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/geom/primitives/amr_faces.h>
#include <zephyr/mesh/euler/eu_face.h>
#include <zephyr/mesh/euler/eu_cell.h>

namespace zephyr::mesh {

EuFace::EuFace(const EuCell &cell, int fid, Direction _dir)
        : m_cell(cell), face_idx(fid), dir(_dir) {

    // Ищем первую определенную грань
    while (face_idx < geom::AmrFaces::max_count && geom().to_skip(dir)) {
        ++face_idx;
    }
}

//AmrFace& EuFace::face() {
//    return m_cell->faces[face_idx];
//}

const AmrFace& EuFace::geom() const {
    return m_cell.geom().faces[face_idx];
}

const Vector3d &EuFace::normal() const {
    return geom().normal;
}

EuCell EuFace::neib() const {
    return geom().is_boundary() ? m_cell : m_cell.neib(face_idx);
}

EuCell EuFace::neighbor() const {
    return geom().is_boundary() ? m_cell : m_cell.neib(face_idx);
}

Boundary EuFace::flag() const {
    return geom().boundary;
}

bool EuFace::is_boundary() const {
    return geom().is_boundary();
}

void EuFace::set_boundary(Boundary flag) {
    AmrFace* face = (AmrFace*)(&geom());
    face->boundary = flag;
}

double EuFace::area() const {
    return geom().area;
}

Vector3d EuFace::center() const {
    return geom().center;
}

EuFaces::EuFaces(const EuCell& cell, Direction _dir)
    : m_cell(cell), dir(_dir) {

}

EuFace EuFaces::begin() const {
    return { m_cell, 0, dir };
}

EuFace EuFaces::end() const {
    return {m_cell, geom::AmrFaces::max_count, dir };
}

} // namespace zephyr::mesh
