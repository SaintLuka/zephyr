#include <zephyr/mesh/storage.h>

#include <zephyr/mesh/face.h>

#include <zephyr/mesh/cell.h>

namespace zephyr { namespace mesh {

IFace::IFace(const ICell &cell, int fid, Direction _dir)
        : m_cell(cell), face_idx(fid), dir(_dir) {

    // Ищем первую определенную грань
    while (face_idx < geom::Faces::max_size && geom().to_skip(dir)) {
        ++face_idx;
    }
}

//Face& IFace::face() {
//    return m_cell->faces[face_idx];
//}

const Face& IFace::geom() const {
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

FaceFlag IFace::flag() const {
    return geom().boundary;
}

bool IFace::is_boundary() const {
    return geom().is_boundary();
}

double IFace::area() const {
    return geom().area;
}

Vector3d IFace::center() const {
    return geom().center(m_cell.geom().vertices, m_cell.geom().dim);
}

IFaces::IFaces(const ICell& cell, Direction _dir)
    : m_cell(cell), dir(_dir) {

}

IFace IFaces::begin() const {
    return { m_cell, 0, dir };
}

IFace IFaces::end() const {
    return { m_cell, geom::Faces::max_size, dir };
}

} // namespace mesh
} // namespace zephyr
