

#include <zephyr/mesh/cell.h>

namespace zephyr { namespace mesh {

Storage::Item safe_iterator(Storage &locals, Storage &aliens, int idx) {
    if (idx < locals.size()) {
        return locals[idx];
    } else {
        idx -= locals.size();
        if (idx < aliens.size()) {
            return aliens[idx];
        } else {
            return locals.end();
        }
    }
}

Storage::Item safe_iterator(Storage &locals, Storage &aliens, const Adjacent &adj) {
    if (adj.ghost < aliens.size()) {
        return aliens[adj.ghost];
    }
    else {
        if (adj.index < locals.size()) {
            return locals[adj.index];
        }
        else {
            return locals.end();
        }
    }
}

ICell::ICell(Storage &_locals, Storage &_aliens, int idx)
    : Storage::Item(safe_iterator(_locals, _aliens, idx)),
      locals(_locals), aliens(_aliens) {

}

ICell::ICell(Storage &_locals, Storage &_aliens, const Adjacent &adj)
    : Storage::Item(safe_iterator(_locals, _aliens, adj)),
      locals(_locals), aliens(_aliens) {

}

ICell ICell::neib(const Face &face) const {
    return { locals, aliens, face.adjacent };
}

ICell ICell::neighbor(const Face &face) const {
    return { locals, aliens, face.adjacent };
}

ICell ICell::neib(const IFace &face) const {
    return { locals, aliens, face.geom().adjacent };
}

ICell ICell::neighbor(const IFace &face) const {
    return { locals, aliens, face.geom().adjacent };
}

ICell ICell::neib(int face_idx) const {
    return { locals, aliens, geom().faces[face_idx].adjacent };
}

ICell ICell::neighbor(int face_idx) const {
    return { locals, aliens, geom().faces[face_idx].adjacent };
}

IFaces ICell::faces() const {
    return { *this };
}

void ICell::print_neibs_info() const {
    print_info();

    std::cout << "\tAll neighbors of cell:\n";
    for (int i = 0; i < geom::Faces::max_size; ++i) {
        auto &face = geom().faces[i];
        if (face.is_undefined() or face.is_boundary()) continue;

        std::cout << "\tNeighbor through the " << side_to_string(geom::Side(i % 6)) << " face (" << i / 6 << "):\n";

        if (face.adjacent.ghost > std::numeric_limits<int>::max()) {
            // Локальная ячейка
            if (face.adjacent.index >= locals.size()) {
                std::cout << "print_cell_info: Wrong connection #1 (It's acceptable for some intermediate refinement stages)";
            }
            else {
                auto neib = locals[face.adjacent.index];
                neib.print_info();
            }
        }
        else {
            // Удаленная ячейка
            if (face.adjacent.ghost >= aliens.size()) {
                std::cout << "print_cell_info: Wrong connection #2 (It's acceptable for some intermediate refinement stages)";
            }
            else {
                auto neib = aliens[face.adjacent.ghost];
                neib.print_info();
            }
        }
    }
}


} // namespace mesh
} // namespace zephyr
