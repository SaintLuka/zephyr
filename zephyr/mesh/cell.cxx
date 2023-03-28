#include <zephyr/mesh/cell.h>

namespace zephyr { namespace mesh {

Storage::iterator safe_iterator(Storage &locals, Storage &aliens, int idx) {
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

Storage::iterator safe_iterator(Storage &locals, Storage &aliens, const Adjacent &adj) {
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
    : iterator(safe_iterator(_locals, _aliens, idx)),
      locals(_locals), aliens(_aliens) {

}

ICell::ICell(Storage &_locals, Storage &_aliens, const Adjacent &adj)
    : iterator(safe_iterator(_locals, _aliens, adj)),
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


} // namespace mesh
} // namespace zephyr
