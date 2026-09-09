#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/geom/primitives/polyhedron.h>

namespace zephyr::mesh {

EuFace_Iter::EuFace_Iter(
        AmrCells *cells, index_t face_idx, index_t face_end,
        AmrCells *ghosts, Direction dir)
        : m_eu_face{cells, face_idx, ghosts},
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
        AmrCells *ghosts,
        Direction dir)
        :
        m_begin(cells,
                cells->faces.offsets[cell_idx],
                cells->faces.offsets[cell_idx + 1],
                ghosts, dir),
        m_end(cells,
              cells->faces.offsets[cell_idx + 1],
              cells->faces.offsets[cell_idx + 1],
              ghosts, dir) { }

geom::Box EuCell::bbox() const {
    return m_cells->bbox(m_index);
}

geom::Polygon EuCell::polygon() const {
    return m_cells->polygon(m_index);
}

geom::Polyhedron EuCell::polyhedron() const {
    return m_cells->polyhedron(m_index);
}

EuCell EuCell::neib(index_t i, index_t j) const {
    z_assert(m_cells->dim() == 2, "neib(i, j) error: not 2D mesh");
    z_assert(m_cells != m_ghosts, "neib(i, j) error: assuming start from local cell");
    z_assert(utils::mpi::single() || m_cells->has_nodes(), "neib(i, j) error: set nodes for distributed mesh");

    bool A = i;
    bool B = j;
    if (A || (!A && B)) {

    }

    // Начинаем с самой ячейки
    AmrCells* cells = m_cells;
    index_t idx = m_index;

    // Сдвинуться на соседа со стороны side
    auto moves = [&, this](Side2D side) {
        z_assert(cells->faces.is_simple(idx, side), "Not structured stencil (Side " + side.to_string() + ")");
        index_t iface = cells->faces.offsets[idx] + side;
        std::tie(cells, idx) = cells->faces.adjacent.get_neib(iface, m_cells, m_ghosts);
    };

    // Переходы направо (ничего не делает при i < 0)
    for (int c = 0; c < i; ++c) {
        moves(Side2D::R);
    }
    // Переходы налево (ничего не делает при i > 0)
    for (int c = 0; c > i; --c) {
        moves(Side2D::L);
    }
    // Переходы вверх (ничего не делает при j < 0)
    for (int c = 0; c < j; ++c) {
        moves(Side2D::T);
    }
    // Переходы вниз (ничего не делает при j > 0)
    for (int c = 0; c > j; --c) {
        moves(Side2D::B);
    }
    return EuCell(cells, idx);
}

EuCell EuCell::neib(index_t i, index_t j, index_t k) const {
    z_assert(m_cells->dim() == 3, "neib(i, j, k) error: not 3D mesh");
    z_assert(m_cells != m_ghosts, "neib(i, j, k) error: assuming start from local cell");
    z_assert(utils::mpi::single() || m_cells->has_nodes(), "neib(i, j, k) error: set nodes for distributed mesh");

    // Начинаем с самой ячейки
    AmrCells* cells = m_cells;
    index_t idx = m_index;

    // Сдвинуться на соседа со стороны side
    auto moves = [&, this](Side3D side) {
        z_assert(cells->faces.is_simple(idx, side), "Not structured stencil (Side " + side.to_string() + ")");
        index_t iface = cells->faces.offsets[idx] + side;
        std::tie(cells, idx) = cells->faces.adjacent.get_neib(iface, m_cells, m_ghosts);
    };

    // Переходы направо (ничего не делает при i < 0)
    for (int c = 0; c < i; ++c) {
        moves(Side3D::R);
    }
    // Переходы налево (ничего не делает при i > 0)
    for (int c = 0; c > i; --c) {
        moves(Side3D::L);
    }
    // Переходы вверх (ничего не делает при j < 0)
    for (int c = 0; c < j; ++c) {
        moves(Side3D::T);
    }
    // Переходы вниз (ничего не делает при j > 0)
    for (int c = 0; c > j; --c) {
        moves(Side3D::B);
    }
    // Переходы вверх (ничего не делает при k < 0)
    for (int c = 0; c < k; ++c) {
        moves(Side3D::F);
    }
    // Переходы вниз (ничего не делает при k > 0)
    for (int c = 0; c > k; --c) {
        moves(Side3D::Z);
    }
    return EuCell(cells, idx);
}

} // namespace zephyr::mesh