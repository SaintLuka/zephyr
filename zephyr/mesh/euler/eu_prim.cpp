#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/geom/primitives/polyhedron.h>

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

geom::Box EuCell::bbox() const {
    return m_cells->bbox(m_index);
}

geom::Polygon EuCell::polygon() const {
    return m_cells->polygon(m_index);
}

geom::Polyhedron EuCell::polyhedron() const {
    return m_cells->polyhedron(m_index);
}

void EuCell::replace(int loc_face) {
    // Переход от alien-ячейки невозможен
    z_assert(m_cells != m_aliens, "Not a local cell #1")
    z_assert(m_cells->rank[m_index] == utils::mpi::rank(), "Not a local cell #2");

    // Индекс правой грани
    index_t iface = m_cells->face_begin[m_index] + loc_face;

    // Массив, в котором находится правая ячейка, индекс ячейки в этом массиве
    std::tie(m_cells, m_index) = m_cells->faces.adjacent.get_neib(iface, m_cells, m_aliens);
}

EuCell EuCell::neib(index_t i, index_t j) const {
    z_assert(m_cells->dim() == 2, "neib(i, j) error: not 2D mesh");

    EuCell neighbor(*this);

    // Переходы направо (ничего не делает при i < 0)
    for (int c = 0; c < i; ++c) {
        z_assert(neighbor.simple_face(Side2D::R), "Not structured stencil (R)");
        neighbor.replace(Side2D::R);
    }
    // Переходы налево (ничего не делает при i > 0)
    for (int c = 0; c > i; --c) {
        z_assert(neighbor.simple_face(Side2D::L), "Not structured stencil (L)");
        neighbor.replace(Side2D::L);
    }
    // Переходы вверх (ничего не делает при j < 0)
    for (int c = 0; c < j; ++c) {
        z_assert(neighbor.simple_face(Side2D::T), "Not structured stencil (T)");
        neighbor.replace(Side2D::T);
    }
    // Переходы вниз (ничего не делает при j > 0)
    for (int c = 0; c > j; --c) {
        z_assert(neighbor.simple_face(Side2D::B), "Not structured stencil (B)");
        neighbor.replace(Side2D::B);
    }
    return neighbor;
}

EuCell EuCell::neib(index_t i, index_t j, index_t k) const {
    z_assert(m_cells->dim() == 3, "neib(i, j, k) error: not 3D mesh");

    EuCell neighbor(*this);

    // Переходы направо (ничего не делает при i < 0)
    for (int c = 0; c < i; ++c) {
        z_assert(neighbor.simple_face(Side3D::R), "Not structured stencil (R)");
        neighbor.replace(Side3D::R);
    }
    // Переходы налево (ничего не делает при i > 0)
    for (int c = 0; c > i; --c) {
        z_assert(neighbor.simple_face(Side3D::L), "Not structured stencil (L)");
        neighbor.replace(Side3D::L);
    }
    // Переходы вверх (ничего не делает при j < 0)
    for (int c = 0; c < j; ++c) {
        z_assert(neighbor.simple_face(Side3D::T), "Not structured stencil (T)");
        neighbor.replace(Side3D::T);
    }
    // Переходы вниз (ничего не делает при j > 0)
    for (int c = 0; c > j; --c) {
        z_assert(neighbor.simple_face(Side3D::B), "Not structured stencil (B)");
        neighbor.replace(Side3D::B);
    }
    // Переходы вверх (ничего не делает при k < 0)
    for (int c = 0; c < k; ++c) {
        z_assert(neighbor.simple_face(Side3D::F), "Not structured stencil (F)");
        neighbor.replace(Side3D::F);
    }
    // Переходы вниз (ничего не делает при k > 0)
    for (int c = 0; c > k; --c) {
        z_assert(neighbor.simple_face(Side3D::Z), "Not structured stencil (Z)");
        neighbor.replace(Side3D::Z);
    }
    return neighbor;
}

} // namespace zephyr::mesh