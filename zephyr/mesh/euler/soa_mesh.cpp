#include <zephyr/mesh/euler/soa_mesh.h>

namespace zephyr::mesh {

void SoaCell::resize(index_t n_cells,
    index_t faces_per_cell, index_t nodes_per_cell) {
    // Поля ячеек по числу ячеек, логично
    rank.resize(n_cells, -1);
    index.resize(n_cells, -1);
    next.resize(n_cells, -1);

    center.resize(n_cells);
    volume.resize(n_cells);
    volume_alt.resize(n_cells);

    // +1 для заключительной
    face_begin.resize(n_cells + 1);
    node_begin.resize(n_cells + 1);

    flag.resize(n_cells);
    b_idx.resize(n_cells);
    z_idx.resize(n_cells);
    level.resize(n_cells);

    // Поля граней и вершин пока так
    faces.resize(n_cells * faces_per_cell);
    verts.resize(n_cells * nodes_per_cell);

    // Поля данных только для ячеек
    data.resize(n_cells);
}

inline double trimin(double x, double y, double z) {
    return std::min(x, std::min(y, z));
}

double QCell::diameter() const {
    if (m_cells->adaptive) {
        if (m_cells->dim == 2) {
            const SqQuad& vertices = m_cells->get_vertices<2>(cell_idx);
            return std::sqrt(std::min(
                    (vertices.vs<+1, 0>() - vertices.vs<-1, 0>()).squaredNorm(),
                    (vertices.vs<0, +1>() - vertices.vs<0, -1>()).squaredNorm()));
        } else {
            const SqCube& vertices = m_cells->get_vertices<3>(cell_idx);
            return std::sqrt(trimin(
                    (vertices.vs<+1, 0, 0>() - vertices.vs<-1, 0, 0>()).squaredNorm(),
                    (vertices.vs<0, +1, 0>() - vertices.vs<0, -1, 0>()).squaredNorm(),
                    (vertices.vs<0, 0, +1>() - vertices.vs<0, 0, -1>()).squaredNorm()));
        }
    }
    else {
        throw std::runtime_error("Not implemented");
    }
}

bool FaceIt::to_skip(Direction dir) const {
    if (m_cells->faces.boundary[face_idx] == Boundary::UNDEFINED) {
        return true;
    }
    switch (dir) {
        case Direction::ANY:
            return false;
        case Direction::X:
            return std::abs(m_cells->faces.normal[face_idx].x()) < 0.7;
        case Direction::Y:
            return std::abs(m_cells->faces.normal[face_idx].y()) < 0.7;
        case Direction::Z:
            return std::abs(m_cells->faces.normal[face_idx].z()) < 0.7;
        default:
            return false;
    }
}

FacesIts QCell::faces(Direction dir) {
    return FacesIts(m_cells,
                    m_cells->face_begin[cell_idx],
                    m_cells->face_begin[cell_idx + 1],
                    dir);
}

bool QFace::is_boundary() const {
    return m_cells->faces.boundary[face_idx] != Boundary::ORDINARY &&
           m_cells->faces.boundary[face_idx] != Boundary::PERIODIC &&
           m_cells->faces.boundary[face_idx] != Boundary::UNDEFINED;
}

bool QFace::is_actual() const {
    return m_cells->faces.boundary[face_idx] != Boundary::UNDEFINED;
}

bool QFace::is_undefined() const {
    return  m_cells->faces.boundary[face_idx] == Boundary::UNDEFINED;
}

void QFace::set_undefined() {
    m_cells->faces.boundary[face_idx] = Boundary::UNDEFINED;
    m_cells->faces.adjacent[face_idx].rank  = -1;
    m_cells->faces.adjacent[face_idx].index = -1;
    m_cells->faces.adjacent[face_idx].alien = -1;
}

Vector3d QFace::symm_point(const Vector3d& p) const {
    return p + 2.0 * (m_cells->faces.center[face_idx] - p).dot(m_cells->faces.normal[face_idx]) * m_cells->faces.normal[face_idx];
}

double QFace::area() const { return m_cells->faces.area[face_idx]; }

double QFace::area(bool axial) const {
    return axial ? m_cells->faces.area_alt[face_idx] : m_cells->faces.area[face_idx];
}

QCell QFace::neib() const {
    return {m_cells, neib_index()};
}

index_t QFace::neib_index() const { return m_cells->faces.adjacent[face_idx].index; }

Vector3d QFace::neib_center() const { return m_cells->center[neib_index()]; }

Boundary QFace::flag() const { return m_cells->faces.boundary[face_idx]; }

const Vector3d& QFace::normal() const { return m_cells->faces.normal[face_idx]; }

const Vector3d& QFace::center() const { return m_cells->faces.center[face_idx]; }

double QFace::x() const { return center().x(); }

inline double QFace::y() const { return center().y(); }

inline double QFace::z() const { return center().z(); }

inline const Adjacent &QFace::adjacent() const {
    return m_cells->faces.adjacent[face_idx];
}


SoaMesh::SoaMesh(EuMesh &mesh)
    : cells(mesh.locals()) {

}

SoaCell::SoaCell(AmrStorage &locals) {
    n_cells = locals.size();

    if (locals.empty()) {
        throw std::runtime_error("SoaMesh: locals is empty");
    }

    dim = locals[0].dim;
    adaptive = locals[0].adaptive;
    linear = locals[0].linear;
    axial = locals[0].axial;

    index_t faces_per_cell = dim < 3 ? 8 : 24;
    index_t nodes_per_cell = dim < 3 ? 9 : 27;

    resize(n_cells, faces_per_cell, nodes_per_cell);


    for (index_t ic = 0; ic < n_cells; ++ic) {
        rank[ic] = locals[ic].rank;
        index[ic] = locals[ic].index;
        next[ic] = locals[ic].next;

        center[ic] = locals[ic].center;
        volume[ic] = locals[ic].volume;
        volume_alt[ic] = locals[ic].volume_alt;
        face_begin[ic] = ic * faces_per_cell;
        node_begin[ic] = ic * nodes_per_cell;

        for (index_t jf = 0; jf < faces_per_cell; ++jf) {
            index_t f_idx = ic * faces_per_cell + jf;

            faces.boundary[f_idx] = locals[ic].faces[jf].boundary;
            faces.adjacent[f_idx] = locals[ic].faces[jf].adjacent;
            faces.normal[f_idx] = locals[ic].faces[jf].normal;
            faces.center[f_idx] = locals[ic].faces[jf].center;
            faces.area[f_idx] = locals[ic].faces[jf].area;
            faces.area_alt[f_idx] = locals[ic].faces[jf].area_alt;
            faces.vertices[f_idx] = locals[ic].faces[jf].vertices;
        }

        for (index_t jn = 0; jn < nodes_per_cell; ++jn) {
            index_t v_idx = ic * nodes_per_cell + jn;
            verts[v_idx] = locals[ic].vertices[jn];
        }
    }
    face_begin[n_cells] = n_cells * faces_per_cell;
    node_begin[n_cells] = n_cells * nodes_per_cell;

    data.resize(n_cells);
}

} // namespace zephyr::mesh