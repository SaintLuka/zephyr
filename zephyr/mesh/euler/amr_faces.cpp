#include <zephyr/geom/cell_type.h>
#include <zephyr/mesh/euler/amr_cells.h>

namespace zephyr::mesh {

using geom::CellType;

void AmrAdjacent::resize(index_t n_faces) {
    rank.resize(n_faces, -1);
    index.resize(n_faces, -1);
    alien.resize(n_faces, -1);
    basic.resize(n_faces, -1);
}

void AmrAdjacent::reserve(index_t n_faces) {
    rank.reserve(n_faces);
    index.reserve(n_faces);
    alien.reserve(n_faces);
    basic.reserve(n_faces);
}

void AmrAdjacent::shrink_to_fit() {
    rank.shrink_to_fit();
    index.shrink_to_fit();
    alien.shrink_to_fit();
    basic.shrink_to_fit();
}

void AmrFaces::resize(index_t n_faces) {
    boundary.resize(n_faces);
    adjacent.resize(n_faces);
    normal.resize(n_faces);
    center.resize(n_faces);
    area.resize(n_faces);
    area_alt.resize(n_faces);
    vertices.resize(n_faces);
}

void AmrFaces::reserve(index_t n_faces) {
    boundary.reserve(n_faces);
    adjacent.reserve(n_faces);
    normal.reserve(n_faces);
    center.reserve(n_faces);
    area.reserve(n_faces);
    area_alt.reserve(n_faces);
    vertices.reserve(n_faces);
}

void AmrFaces::shrink_to_fit() {
    boundary.shrink_to_fit();
    adjacent.shrink_to_fit();
    normal.shrink_to_fit();
    center.shrink_to_fit();
    area.shrink_to_fit();
    area_alt.shrink_to_fit();
    vertices.shrink_to_fit();
}

void AmrFaces::insert(index_t iface, CellType ctype, int count) {
    int n_faces = -1;
    switch (ctype) {
        case CellType::AMR2D:
            vertices[iface + Side2D::L] = Side2D::L.sf();
            vertices[iface + Side2D::R] = Side2D::R.sf();
            vertices[iface + Side2D::B] = Side2D::B.sf();
            vertices[iface + Side2D::T] = Side2D::T.sf();
            n_faces = 8;
            break;

        case CellType::AMR3D:
            vertices[iface + Side3D::L] = Side3D::L.sf();
            vertices[iface + Side3D::R] = Side3D::R.sf();
            vertices[iface + Side3D::B] = Side3D::B.sf();
            vertices[iface + Side3D::T] = Side3D::T.sf();
            vertices[iface + Side3D::X] = Side3D::X.sf();
            vertices[iface + Side3D::F] = Side3D::F.sf();
            n_faces = 24;
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
                vertices[iface + i] = {i, (i + 1) % count, -1, -1};
            }
            n_faces = count;
            break;

        case CellType::POLYHEDRON:
            break;

        default:
            std::string message = "BFaces::BFaces error(): not implemented for current CellType";
            std::cerr << message << "\n";
            throw std::runtime_error(message);
    }

    for (int i = 0; i < n_faces; ++i) {
        set_undefined(iface + i);
    }
}

bool AmrFaces::to_skip(index_t iface, Direction dir) const {
    if (boundary[iface] == Boundary::UNDEFINED) {
        return true;
    }
    switch (dir) {
        case Direction::ANY: return false;
        case Direction::X: return std::abs(normal[iface].x()) < 0.7;
        case Direction::Y: return std::abs(normal[iface].y()) < 0.7;
        case Direction::Z: return std::abs(normal[iface].z()) < 0.7;
        default: return false;
    }
}

} // namespace zephyr::mesh