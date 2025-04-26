#include <zephyr/mesh/primitives/Side3D.h>
#include <zephyr/mesh/primitives/amr_cell.h>
#include <zephyr/mesh/euler/eu_cell.h>

namespace zephyr::mesh {

AmrStorage::Iterator safe_iterator(AmrStorage &locals, AmrStorage &aliens, int idx) {
    if (idx <= locals.size()) {
        return locals.iterator(idx);
    } else {
        idx -= locals.size();
        if (idx < aliens.size()) {
            return aliens.iterator(idx);
        } else {
            return locals.end();
        }
    }
}

AmrStorage::Iterator safe_iterator(AmrStorage &locals, AmrStorage &aliens, const Adjacent &adj) {
    if (adj.alien >= 0) {
        return aliens.iterator(adj.alien);
    }
    else {
        if (adj.index < locals.size()) {
            return locals.iterator(adj.index);
        }
        else {
            return locals.end();
        }
    }
}

EuCell::EuCell(AmrStorage &_locals, AmrStorage &_aliens, size_t idx)
    : m_it(safe_iterator(_locals, _aliens, idx)),
      m_locals(_locals), m_aliens(_aliens) {

}

EuCell::EuCell(AmrStorage &_locals, AmrStorage &_aliens, const Adjacent &adj)
    : m_it(safe_iterator(_locals, _aliens, adj)),
      m_locals(_locals), m_aliens(_aliens) {

}

EuCell EuCell::locals(size_t idx) {
    return { m_locals, m_aliens, idx };
}

EuCell EuCell::neib(const BFace &face) const {
    return { m_locals, m_aliens, face.adjacent };
}

EuCell EuCell::neib(const EuFace &face) const {
    return { m_locals, m_aliens, face.adjacent() };
}

EuCell EuCell::neib(int face_idx) const {
    return { m_locals, m_aliens, m_it->faces[face_idx].adjacent };
}

const Byte* EuCell::neib_data(const BFace& face) const {
    if (face.is_undefined() || face.is_boundary()) {
        return m_it->data();
    }
    const auto &adj = face.adjacent;

    if (adj.alien >= 0) {
        return m_aliens.data_at(adj.alien);
    } else if (adj.index < m_locals.size()) {
        return m_locals.data_at(adj.index);
    } else {
        return m_it->data();
    }
}

EuFace EuCell::face(Side3D side) const {
    return EuFace(*this, side);
}

EuFaces EuCell::faces(Direction dir) const {
    return { *this, dir };
}

void EuCell::print_neibs_info() const {
    m_it->print_info();

    std::cout << "\tAll neighbors of cell:\n";
    for (int i = 0; i < BFaces::max_count; ++i) {
        auto &face = m_it->faces[i];
        if (face.is_undefined() or face.is_boundary()) continue;

        std::cout << "\tNeighbor through the " << side_to_string(i) << ":\n";

        if (face.adjacent.alien > std::numeric_limits<int>::max()) {
            // Локальная ячейка
            if (face.adjacent.index >= m_locals.size()) {
                std::cout << "print_cell_info: Wrong connection #1 (It's acceptable for some intermediate refinement stages)";
            }
            else {
                auto& neib = m_locals[face.adjacent.index];
                neib.print_info();
            }
        }
        else {
            // Удаленная ячейка
            if (face.adjacent.alien >= m_aliens.size()) {
                std::cout << "print_cell_info: Wrong connection #2 (It's acceptable for some intermediate refinement stages)";
            }
            else {
                auto& neib = m_aliens[face.adjacent.alien];
                neib.print_info();
            }
        }
    }
}

} // namespace zephyr::mesh
