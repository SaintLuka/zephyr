#include <cstring>

#include <zephyr/geom/primitives/amr_cell.h>
#include <zephyr/geom/primitives/amr_faces.h>
#include <zephyr/mesh/face.h>
#include <zephyr/mesh/cell.h>

namespace zephyr::mesh {

IFace::IFace(const ICell &cell, int fid, Direction _dir)
        : m_cell(cell), face_idx(fid), dir(_dir) {

    // Ищем первую определенную грань
    while (face_idx < geom::AmrFaces::max_count && geom().to_skip(dir)) {
        ++face_idx;
    }
}

//AmrFace& IFace::face() {
//    return m_cell->faces[face_idx];
//}

const AmrFace& IFace::geom() const {
    return m_cell.geom().faces[face_idx];
}

const Vector3d &IFace::normal() const {
    return geom().normal;
}

ICell IFace::neib() const {
    return geom().is_boundary() ? m_cell : m_cell.neib(face_idx);
}

ICell IFace::neighbor() const {
    return geom().is_boundary() ? m_cell : m_cell.neib(face_idx);
}

Boundary IFace::flag() const {
    return geom().boundary;
}

bool IFace::is_boundary() const {
    return geom().is_boundary();
}

void IFace::set_boundary(Boundary flag) {
    AmrFace* face = (AmrFace*)(&geom());
    face->boundary = flag;
}

double IFace::area() const {
    return geom().area;
}

Vector3d IFace::center() const {
    return geom().center;
}

IFaces::IFaces(const ICell& cell, Direction _dir)
    : m_cell(cell), dir(_dir) {

}

IFace IFaces::begin() const {
    return { m_cell, 0, dir };
}

IFace IFaces::end() const {
    return {m_cell, geom::AmrFaces::max_count, dir };
}

} // namespace zephyr::mesh
