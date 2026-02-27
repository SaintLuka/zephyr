#include <zephyr/geom/cell_type.h>
#include <zephyr/mesh/euler/amr_cells.h>

#include "zephyr/geom/generator/array2d.h"

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
            z_assert(count == 8 || count < 0, "AmrFaces::insert: bad size 2D");
            vertices[iface + Side2D::L] = Side2D::L.sf();
            vertices[iface + Side2D::R] = Side2D::R.sf();
            vertices[iface + Side2D::B] = Side2D::B.sf();
            vertices[iface + Side2D::T] = Side2D::T.sf();
            n_faces = 8;
            break;

        case CellType::AMR3D:
            z_assert(count == 24 || count < 0, "AmrFaces::insert: bad size 3D");
            vertices[iface + Side3D::L] = Side3D::L.sf();
            vertices[iface + Side3D::R] = Side3D::R.sf();
            vertices[iface + Side3D::B] = Side3D::B.sf();
            vertices[iface + Side3D::T] = Side3D::T.sf();
            vertices[iface + Side3D::Z] = Side3D::Z.sf();
            vertices[iface + Side3D::F] = Side3D::F.sf();
            n_faces = 24;
            break;

        case CellType::TRIANGLE:
            z_assert(count == 3 || count < 0, "AmrFaces::insert: bad size TRIANGLE");
            vertices[iface+0] = {0, 1, -1, -1, -1, -1, -1, -1};
            vertices[iface+1] = {1, 2, -1, -1, -1, -1, -1, -1};
            vertices[iface+2] = {2, 0, -1, -1, -1, -1, -1, -1};       
            n_faces = 3;
            break;
            
        case CellType::QUAD:
            // Необычный порядок граней
            z_assert(count == 4 || count < 0, "AmrFaces::insert: bad size QUAD");
            vertices[iface+0] = {0, 3, -1, -1, -1, -1, -1, -1};
            vertices[iface+1] = {1, 2, -1, -1, -1, -1, -1, -1};
            vertices[iface+2] = {0, 1, -1, -1, -1, -1, -1, -1};
            vertices[iface+3] = {3, 2, -1, -1, -1, -1, -1, -1};
            n_faces = 4;
            break;

        case CellType::POLYGON:
            if (count < 0) {
                throw std::runtime_error("AmrFaces::insert error: set argument 'count' with CellType::POLYGON");
            }
            for (int i = 0; i < count; ++i) {
                vertices[iface + i].fill(-1);
                vertices[iface + i][0] = static_cast<short>(i);
                vertices[iface + i][1] = static_cast<short>((i + 1) % count);
            }
            n_faces = count;
            break;

        case CellType::HEXAHEDRON:
            z_assert(count == 6 || count < 0, "AmrFaces::insert: bad size HEXAHEDRON");
            vertices[iface + Side3D::L] = {7, 3, 0, 4, -1, -1, -1, -1};
            vertices[iface + Side3D::R] = {1, 2, 6, 5, -1, -1, -1, -1};
            vertices[iface + Side3D::B] = {0, 1, 5, 4, -1, -1, -1, -1};
            vertices[iface + Side3D::T] = {2, 3, 7, 6, -1, -1, -1, -1};
            vertices[iface + Side3D::B] = {3, 2, 1, 0, -1, -1, -1, -1};
            vertices[iface + Side3D::F] = {4, 5, 6, 7, -1, -1, -1, -1};
            n_faces = 6;
            break;

        case CellType::POLYHEDRON:
        default:
            throw std::runtime_error("AmrFaces::insert: not implemented for current CellType");
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

template <class T>
void reorder(std::vector<T>& field, index_t iface) {
    // 0 -> 2, 2 -> 3, 3 -> 0
    T f0 = field[iface];
    field[iface] = field[iface + 3];
    field[iface + 3] = field[iface + 2];
    field[iface + 2] = f0;
}

void AmrFaces::reorder_quad_faces(index_t iface) {
    reorder(adjacent.rank, iface);
    reorder(adjacent.index, iface);
    reorder(adjacent.alien, iface);

    reorder(boundary, iface);
    reorder(normal, iface);
    reorder(center, iface);
    reorder(area, iface);
    reorder(vertices, iface);

    if (!area_alt.empty()) {
        reorder(area_alt, iface);
    }
}

} // namespace zephyr::mesh