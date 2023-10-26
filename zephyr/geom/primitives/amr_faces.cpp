#include <zephyr/geom/primitives/amr_faces.h>

namespace zephyr::geom {

AmrFaces::AmrFaces(int dim) {
    if (dim < 3) {
        m_list[Side::L].vertices = face_indices::sf<2, Side::L>();
        m_list[Side::R].vertices = face_indices::sf<2, Side::R>();
        m_list[Side::B].vertices = face_indices::sf<2, Side::B>();
        m_list[Side::T].vertices = face_indices::sf<2, Side::T>();
    } else {
        m_list[Side::L].vertices = face_indices::sf<3, Side::L>();
        m_list[Side::R].vertices = face_indices::sf<3, Side::R>();
        m_list[Side::B].vertices = face_indices::sf<3, Side::B>();
        m_list[Side::T].vertices = face_indices::sf<3, Side::T>();
        m_list[Side::X].vertices = face_indices::sf<3, Side::X>();
        m_list[Side::F].vertices = face_indices::sf<3, Side::F>();
    }
}

int AmrFaces::count() const {
    int count = 0;
    for (auto& face: m_list) {
        if (face.is_actual()) {
            ++count;
        }
    }
    return count;
}

} // namespace zephyr::geom