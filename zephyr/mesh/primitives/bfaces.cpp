#include <zephyr/mesh/primitives/bfaces.h>

namespace zephyr::mesh {

BFaces::BFaces(CellType ctype, int count) {
    switch (ctype) {
        case CellType::AMR2D:
            m_faces[Side3D::L].vertices = Side2D::L.sf();
            m_faces[Side3D::R].vertices = Side2D::R.sf();
            m_faces[Side3D::B].vertices = Side2D::B.sf();
            m_faces[Side3D::T].vertices = Side2D::T.sf();
            break;

        case CellType::AMR3D:
            m_faces[Side3D::L].vertices = Side3D::L.sf();
            m_faces[Side3D::R].vertices = Side3D::R.sf();
            m_faces[Side3D::B].vertices = Side3D::B.sf();
            m_faces[Side3D::T].vertices = Side3D::T.sf();
            m_faces[Side3D::X].vertices = Side3D::X.sf();
            m_faces[Side3D::F].vertices = Side3D::F.sf();
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
                m_faces[i].vertices = {i, (i + 1) % count, -1, -1};
            }
            break;

        case CellType::POLYHEDRON:
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

} // namespace zephyr::mesh