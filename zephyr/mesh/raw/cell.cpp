#include <zephyr/mesh/cell.h>
#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/geom/primitives/polyhedron.h>

namespace zephyr::mesh {

Face_Iter::Face_Iter(
        RawCells *cells, index_t face_idx, index_t face_end,
        RawCells *ghosts, Direction dir)
        : face_{cells, face_idx, ghosts},
          face_end_(face_end),
          dir_(dir) {

    while (face_.face_idx_ < face_end_ && to_skip(dir_)) {
        face_.face_idx_ += 1;
    }
}

Face_Iter &Face_Iter::operator++() {
    do {
        face_.face_idx_ += 1;
    } while (face_.face_idx_ < face_end_ && to_skip(dir_));
    return *this;
}

bool Face_Iter::operator!=(const Face_Iter &face) const {
    return face_.face_idx_ != face.face_.face_idx_;
}

bool Face_Iter::to_skip(Direction dir) const {
    return face_.cells_->faces.to_skip(face_.face_idx_, dir);
}

FacesRange::FacesRange(
        RawCells *cells,
        index_t cell_idx,
        RawCells *ghosts,
        Direction dir)
        :
        begin_(cells,
                cells->faces.offsets[cell_idx],
                cells->faces.offsets[cell_idx + 1],
                ghosts, dir),
        end_(cells,
              cells->faces.offsets[cell_idx + 1],
              cells->faces.offsets[cell_idx + 1],
              ghosts, dir) { }

geom::Box Cell::bbox() const {
    return cells_->bbox(index_);
}

geom::Polygon Cell::polygon() const {
    return cells_->polygon(index_);
}

geom::Polyhedron Cell::polyhedron() const {
    return cells_->polyhedron(index_);
}

Cell Cell::neib(index_t i, index_t j) const {
    z_assert(cells_->dim() == 2, "neib(i, j) error: not 2D mesh");
    z_assert(cells_ != ghosts_, "neib(i, j) error: assuming start from local cell");
    z_assert(utils::mpi::single() || cells_->has_nodes(), "neib(i, j) error: set nodes for distributed mesh");

    // Начинаем с самой ячейки
    RawCells* cells = cells_;
    index_t idx = index_;

    // Сдвинуться на соседа со стороны side
    auto moves = [&, this](Side2D side) {
        z_assert(cells->faces.is_simple(idx, side), "Not structured stencil (Side " + side.to_string() + ")");
        index_t iface = cells->faces.offsets[idx] + side;
        std::tie(cells, idx) = cells->faces.adjacent.get_neib(iface, cells_, ghosts_);
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
    return Cell(cells, idx);
}

Cell Cell::neib(index_t i, index_t j, index_t k) const {
    z_assert(cells_->dim() == 3, "neib(i, j, k) error: not 3D mesh");
    z_assert(cells_ != ghosts_, "neib(i, j, k) error: assuming start from local cell");
    z_assert(utils::mpi::single() || cells_->has_nodes(), "neib(i, j, k) error: set nodes for distributed mesh");

    // Начинаем с самой ячейки
    RawCells* cells = cells_;
    index_t idx = index_;

    // Сдвинуться на соседа со стороны side
    auto moves = [&, this](Side3D side) {
        z_assert(cells->faces.is_simple(idx, side), "Not structured stencil (Side " + side.to_string() + ")");
        index_t iface = cells->faces.offsets[idx] + side;
        std::tie(cells, idx) = cells->faces.adjacent.get_neib(iface, cells_, ghosts_);
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
    return Cell(cells, idx);
}

bool Cell::local_neibs() const {
    for (index_t iface: cells_->faces.range(index_)) {
        auto flag = cells_->faces.boundary[iface];
        if (flag == geom::Boundary::INNER || flag == geom::Boundary::PERIODIC) {
            if (cells_->faces.adjacent.is_ghost(iface)) {
                return false;
            }
        }
    }
    return true;
}

// Сосед того же уровня по стороне side
template<int dim>
bool good_neibs(const Cell& cell, Side<dim> side, int lvl) {
    return cell.simple_face(side) && cell.face(side).neib_level() == lvl;
}

// Локальный сосед того же уровня по стороне side
template<int dim>
bool good_local_neibs(const Cell& cell, Side<dim> side, int lvl) {
    return cell.simple_face(side) && cell.face(side).local_neib() && cell.face(side).neib_level() == lvl;
}

bool Cell::extended_stencil() const {
    int lvl = cells_->level[index_];

    if (dim() < 3) {
        if (!good_local_neibs(*this, Side2D::L, lvl)) return false;
        if (!good_local_neibs(*this, Side2D::R, lvl)) return false;
        if (!good_local_neibs(*this, Side2D::B, lvl)) return false;
        if (!good_local_neibs(*this, Side2D::T, lvl)) return false;

        auto neib_L = face(Side2D::L).neib();
        if (!good_neibs(neib_L, Side2D::B, lvl)) return false;
        if (!good_neibs(neib_L, Side2D::T, lvl)) return false;

        auto neib_R = face(Side2D::R).neib();
        if (!good_neibs(neib_R, Side2D::B, lvl)) return false;
        if (!good_neibs(neib_R, Side2D::T, lvl)) return false;

        return true;
    }
    else {
        // Хорошие локальные соседи через 6 граней
        for (auto side: Side3D::items()) {
            if (!good_local_neibs(*this, side, lvl)) return false;
        }

        // Хорошие локальные соседи второго уровня
        for (auto side_x: {Side3D::L, Side3D::R}) {
            auto neib_x = face(side_x).neib();
            for (auto side: {Side3D::B, Side3D::T, Side3D::Z, Side3D::F}) {
                if (!good_local_neibs(neib_x, side, lvl)) return false;
            }
        }

        for (auto side_y: {Side3D::B, Side3D::T}) {
            auto neib_y = face(side_y).neib();
            for (auto side: {Side3D::L, Side3D::R, Side3D::Z, Side3D::F}) {
                if (!good_local_neibs(neib_y, side, lvl)) return false;
            }
        }

        for (auto side_z: {Side3D::Z, Side3D::F}) {
            auto neib_z = face(side_z).neib();
            for (auto side: {Side3D::L, Side3D::R, Side3D::B, Side3D::T}) {
                if (!good_local_neibs(neib_z, side, lvl)) return false;
            }
        }

        // Хорошие угловые соседи
        for (auto side: {Side3D::Z, Side3D::F}) {
            if (!good_neibs(neib(-1, -1, 0), side, lvl)) return false;
            if (!good_neibs(neib(-1, +1, 0), side, lvl)) return false;
            if (!good_neibs(neib(+1, -1, 0), side, lvl)) return false;
            if (!good_neibs(neib(+1, +1, 0), side, lvl)) return false;
        }

        return true;
    }
    return true;
}

} // namespace zephyr::mesh