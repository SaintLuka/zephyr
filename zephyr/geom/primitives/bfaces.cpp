#include <zephyr/geom/primitives/bfaces.h>

namespace zephyr::geom {

BFaces::BFaces(CellType ctype, int count) {
    switch (ctype) {
        case CellType::AMR2D:
            m_faces[Side::L].vertices = face_indices::sf<2, Side::L>();
            m_faces[Side::R].vertices = face_indices::sf<2, Side::R>();
            m_faces[Side::B].vertices = face_indices::sf<2, Side::B>();
            m_faces[Side::T].vertices = face_indices::sf<2, Side::T>();
            break;

        case CellType::AMR3D:
            m_faces[Side::L].vertices = face_indices::sf<3, Side::L>();
            m_faces[Side::R].vertices = face_indices::sf<3, Side::R>();
            m_faces[Side::B].vertices = face_indices::sf<3, Side::B>();
            m_faces[Side::T].vertices = face_indices::sf<3, Side::T>();
            m_faces[Side::X].vertices = face_indices::sf<3, Side::X>();
            m_faces[Side::F].vertices = face_indices::sf<3, Side::F>();
            break;

        case CellType::TRIANGLE:
        case CellType::QUAD:
            // определяем count и идем на CellType::POLYGON
            count = ctype == CellType::TRIANGLE ? 3 : 4;

        case CellType::POLYGON:
            if (count < 0) {
                std::string message = "BFaces::BFaces error(): set argument 'count' with CellType::POLYGON";
                std::cerr << message << "\n";
                throw std::runtime_error(message);
            }
            for (int i = 0; i < count; ++i) {
                m_faces[i].vertices = {i, (i + 1) % count,
                                       face_indices::undef(),
                                       face_indices::undef()};
            }
            break;

        default:
            std::string message = "BFaces::BFaces error(): not implemented for current CellType";
            std::cerr << message << "\n";
            throw std::runtime_error(message);
    }
}

int BFaces::count() const {
    int count = 0;
    for (auto& face: m_faces) {
        if (face.is_actual()) {
            ++count;
        }
    }
    return count;
}

} // namespace zephyr::geom