#include <zephyr/geom/primitives/side.h>
#include <zephyr/geom/primitives/amr_cell.h>

namespace zephyr::geom {

AmrCell::AmrCell(const Quad& quad)
    : Element(0, 0), dim(2),
    adaptive(true), linear(true),
    vertices(quad), faces(CellType::AMR2D),
    b_idx(-1), z_idx(-1), level(0), flag(0) {

    double area = quad.area();

    size = std::sqrt(area);
    center = quad.centroid(area);

    for (int i = 0; i < 4; ++i) {
        Line vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]]
        };

        faces[i].area     = vs.length();
        faces[i].center   = vs.center();
        faces[i].normal   = vs.normal(center);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

AmrCell::AmrCell(const SqQuad& quad)
    : AmrCell(quad.reduce()) {
    //std::cerr << "Nonlinear AmrCells is not supported\n";
}

AmrCell::AmrCell(const Cube& cube)
    : Element(0, 0), dim(3),
    adaptive(true), linear(true),
    vertices(cube), faces(CellType::AMR3D),
    b_idx(-1), z_idx(-1), level(0), flag(0) {

    double volume = cube.volume();
    size = std::cbrt(volume);
    center = cube.centroid(volume);

    for (int i = 0; i < 6; ++i) {
        Quad vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]],
                vertices[faces[i].vertices[2]],
                vertices[faces[i].vertices[3]]
        };

        faces[i].area     = vs.area();
        faces[i].center   = vs.center();
        faces[i].normal   = vs.normal(center);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

AmrCell::AmrCell(const SqCube& cube)
    : AmrCell(cube.reduce()) {
    //std::cerr << "Nonlinear AmrCells is not supported\n";
}

AmrCell::AmrCell(const Polygon& poly)
        : Element(0, 0), dim(2),
          adaptive(false), linear(true),
          vertices(poly), faces(CellType::POLYGON, poly.size()),
          b_idx(-1), z_idx(-1), level(0), flag(0) {

    Vector3d C = poly.center();
    double area = poly.area(C);

    size = std::sqrt(area);
    center = poly.centroid(area);

    for (int i = 0; i < poly.size(); ++i) {
        Line vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]]
        };

        faces[i].area     = vs.length();
        faces[i].center   = vs.center();
        faces[i].normal   = vs.normal(center);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

double AmrCell::volume() const {
    return size * (dim < 3 ? size : size * size);
}

PolygonS<8> AmrCell::polygon() const {
    if (dim > 2) {
        throw std::runtime_error("AmrCell::polygon() error #1");
    }

    PolygonS<8> poly;
    int n_nodes = 0;

    if (adaptive) {
        poly[n_nodes++] = vertices.vs<-1, -1>();
        if (!linear && faces[Side::BOTTOM1].is_actual()) {
            poly[n_nodes++] = vertices.vs<0, -1>();
        }
        poly[n_nodes++] = vertices.vs<+1, -1>();
        if (!linear && faces[Side::RIGHT1].is_actual()) {
            poly[n_nodes++] = vertices.vs<+1, 0>();
        }
        poly[n_nodes++] = vertices.vs<+1, +1>();
        if (!linear && faces[Side::TOP1].is_actual()) {
            poly[n_nodes++] = vertices.vs<0, +1>();
        }
        poly[n_nodes++] = vertices.vs<-1, +1>();
        if (!linear && faces[Side::LEFT1].is_actual()) {
            poly[n_nodes++] = vertices.vs<-1, 0>();
        }

    } else {
        while (!vertices[n_nodes].hasNaN() && n_nodes < 8) {
            poly[n_nodes] = vertices[n_nodes];
            ++n_nodes;
        }
    }

    poly.resize(n_nodes);
    return poly;
}

void AmrCell::mark_actual_nodes(int mark) {
    static_assert(BNodes::max_count == BVertices::max_count);

    nodes.clear();
    for (const BFace &face: faces) {
        for (int k = 0; k < BFace::max_size; ++k) {
            int iv = face.vertices[k];
            if (iv >= 0) {
                nodes[iv] = -13;
            }
            else {
                break;
            }
        }
    }
}

} // namespace zephyr::geom