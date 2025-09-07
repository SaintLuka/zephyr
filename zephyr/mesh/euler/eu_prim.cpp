#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/geom/primitives/polygon.h>

namespace zephyr::mesh {

EuFace_Iter::EuFace_Iter(
        AmrCells *cells, index_t face_idx, index_t face_end,
        AmrCells *aliens, Direction dir)
        : m_eu_face{cells, face_idx, aliens},
          m_face_end(face_end),
          m_dir(dir) {

    while (m_eu_face.m_face_idx < m_face_end && to_skip(m_dir)) {
        m_eu_face.m_face_idx += 1;
    }
}

EuFace_Iter &EuFace_Iter::operator++() {
    do {
        m_eu_face.m_face_idx += 1;
    } while (m_eu_face.m_face_idx < m_face_end && to_skip(m_dir));
    return *this;
}

bool EuFace_Iter::operator!=(const EuFace_Iter &face) const {
    return m_eu_face.m_face_idx != face.m_eu_face.m_face_idx;
}

bool EuFace_Iter::to_skip(Direction dir) const {
    return m_eu_face.m_cells->faces.to_skip(m_eu_face.m_face_idx, dir);
}

EuFaces::EuFaces(
        AmrCells *cells,
        index_t cell_idx,
        AmrCells *aliens,
        Direction dir)
        :
        m_begin(cells,
                cells->face_begin[cell_idx],
                cells->face_begin[cell_idx + 1],
                aliens, dir),
        m_end(cells,
              cells->face_begin[cell_idx + 1],
              cells->face_begin[cell_idx + 1],
              aliens, dir) { }

geom::Polygon EuCell::polygon() const {
    return m_cells->polygon(m_index);
}

} // namespace zephyr::mesh