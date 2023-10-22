#include <zephyr/geom/primitives/amr_cell.h>

namespace zephyr::geom {

void AmrCell::setup_vertices(const Quad& quad) {
    vertices.set_undefined();

    // Базовые вершины
    vertices[iv(0, 0)] = quad.vs<-1, -1>();
    vertices[iv(0, 1)] = quad.vs<-1, +1>();
    vertices[iv(1, 0)] = quad.vs<+1, -1>();
    vertices[iv(1, 1)] = quad.vs<+1, +1>();

    // Вершины на гранях
    vertices[iww(0, 1)] = quad(-1.0, 0.0);
    vertices[iww(2, 1)] = quad(+1.0, 0.0);
    vertices[iww(1, 0)] = quad(0.0, -1.0);
    vertices[iww(1, 2)] = quad(0.0, +1.0);

    // Вершина в центре
    vertices[iww(1, 1)] = quad(0.0, 0.0);
}

void AmrCell::setup_vertices(const SqQuad& vlist) {
    vertices.set_undefined();

    // Базовые вершины
    vertices[iww(0, 0)] = vlist.vs<-1, -1>();
    vertices[iww(0, 2)] = vlist.vs<-1, +1>();
    vertices[iww(2, 0)] = vlist.vs<+1, -1>();
    vertices[iww(2, 2)] = vlist.vs<+1, +1>();

    // Вершины на гранях
    vertices[iww(0, 1)] = vlist.vs<-1, 0>();
    vertices[iww(2, 1)] = vlist.vs<+1, 0>();
    vertices[iww(1, 0)] = vlist.vs<0, -1>();
    vertices[iww(1, 2)] = vlist.vs<0, +1>();

    // Вершина в центре
    vertices[iww(1, 1)] = vlist.vs<0, 0>();
}

void AmrCell::setup_vertices(const Cube& cube) {
    vertices.set_undefined();

    // Скопировать базовые вершины
    vertices[iv(0, 0, 0)] = cube.vs<-1, -1, -1>();
    vertices[iv(1, 0, 0)] = cube.vs<+1, -1, -1>();
    vertices[iv(0, 1, 0)] = cube.vs<-1, +1, -1>();
    vertices[iv(1, 1, 0)] = cube.vs<+1, +1, -1>();
    vertices[iv(0, 0, 1)] = cube.vs<-1, -1, +1>();
    vertices[iv(1, 0, 1)] = cube.vs<+1, -1, +1>();
    vertices[iv(0, 1, 1)] = cube.vs<-1, +1, +1>();
    vertices[iv(1, 1, 1)] = cube.vs<+1, +1, +1>();

    // Вершины на ребрах
    vertices[iww(0, 1, 0)] = cube(-1.0,  0.0, -1.0);
    vertices[iww(2, 1, 0)] = cube(+1.0,  0.0, -1.0);
    vertices[iww(1, 0, 0)] = cube( 0.0, -1.0, -1.0);
    vertices[iww(1, 2, 0)] = cube( 0.0, +1.0, -1.0);
    vertices[iww(0, 0, 1)] = cube(-1.0, -1.0,  0.0);
    vertices[iww(2, 0, 1)] = cube(+1.0, -1.0,  0.0);
    vertices[iww(0, 2, 1)] = cube(-1.0, +1.0,  0.0);
    vertices[iww(2, 2, 1)] = cube(+1.0, +1.0,  0.0);
    vertices[iww(0, 1, 2)] = cube(-1.0,  0.0, +1.0);
    vertices[iww(2, 1, 2)] = cube(+1.0,  0.0, +1.0);
    vertices[iww(1, 0, 2)] = cube( 0.0, -1.0, +1.0);
    vertices[iww(1, 2, 2)] = cube( 0.0, +1.0, +1.0);

    // Вершны на гранях
    vertices[iww(0, 1, 1)] = cube(-1.0,  0.0,  0.0);
    vertices[iww(2, 1, 1)] = cube(+1.0,  0.0,  0.0);
    vertices[iww(1, 0, 1)] = cube( 0.0, -1.0,  0.0);
    vertices[iww(1, 2, 1)] = cube( 0.0, +1.0,  0.0);
    vertices[iww(1, 1, 0)] = cube( 0.0,  0.0, -1.0);
    vertices[iww(1, 1, 2)] = cube( 0.0,  0.0, +1.0);

    // Вершина в центре
    vertices[iww(1, 1, 1)] = cube( 0.0,  0.0,  0.0);
}

void AmrCell::build2D(const Quad& verts) {
    dim = 2;

    double area = verts.area();
    Vector3d centroid = verts.centroid(area);

    coords = centroid;
    size = std::sqrt(area);

    setup_vertices(verts);

    faces.set_undefined();

    faces[Side::L].vertices = face_indices<2, Side::L>();
    faces[Side::R].vertices = face_indices<2, Side::R>();
    faces[Side::B].vertices = face_indices<2, Side::B>();
    faces[Side::T].vertices = face_indices<2, Side::T>();

    for (int i = 0; i < 4; ++i) {
        Line vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]]
        };

        faces[i].area     = vs.length();
        faces[i].normal   = vs.normal(centroid);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

void AmrCell::build2D(const SqQuad& verts) {
    dim = 2;

    double area = verts.area();
    Vector3d centroid = verts.centroid(area);

    coords = centroid;
    size = std::sqrt(area);

    setup_vertices(verts);

    faces.set_undefined();

    faces[Side::L].vertices = face_indices<2, Side::L>();
    faces[Side::R].vertices = face_indices<2, Side::R>();
    faces[Side::B].vertices = face_indices<2, Side::B>();
    faces[Side::T].vertices = face_indices<2, Side::T>();

    for (int i = 0; i < 4; ++i) {
        Line vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]]
        };

        faces[i].area     = vs.length();
        faces[i].normal   = vs.normal(centroid);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

void AmrCell::build2D(const VerticesList& verts) {
    dim = 2;

    double area = geom::area(verts);
    Vector3d centroid = geom::centroid(verts, area);

    coords = centroid;
    size = std::sqrt(area);

    faces.set_undefined();
    vertices.set_undefined();

    // Используем базисные точки
    int n_points = verts.size();
    for (int i = 0; i < n_points; ++i) {
        vertices[i] = verts[i];
    }

    short undef = undef_index();

    for (int i = 0; i < n_points; ++i) {
        int j = (i + 1) % n_points;

        faces[i].vertices = {short(i), short(j), undef, undef};

        Line vs = {verts[i], verts[j]};

        faces[i].area     = vs.length();
        faces[i].normal   = vs.normal(centroid);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

void AmrCell::build3D(const Cube& verts) {
    dim = 3;

    double volume = verts.volume();
    size = std::cbrt(volume);

    Vector3d centroid = verts.centroid(volume);
    coords = centroid;

    setup_vertices(verts);

    faces.set_undefined();

    faces[Side::L].vertices = face_indices<3, Side::L>();
    faces[Side::R].vertices = face_indices<3, Side::R>();
    faces[Side::B].vertices = face_indices<3, Side::B>();
    faces[Side::T].vertices = face_indices<3, Side::T>();
    faces[Side::X].vertices = face_indices<3, Side::X>();
    faces[Side::F].vertices = face_indices<3, Side::F>();

    for (int i = 0; i < 6; ++i) {
        Quad vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]],
                vertices[faces[i].vertices[2]],
                vertices[faces[i].vertices[3]]
        };

        faces[i].area     = vs.area();
        faces[i].normal   = vs.normal(centroid);
        faces[i].boundary = Boundary::ORDINARY;
    }
}

AmrCell::AmrCell() :
    Element(0, 0),
    b_idx(0),
    z_idx(0),
    next(0),
    level(0),
    flag(0),
    dim(0),
    size(0.0) {

}

AmrCell::AmrCell(const Quad& verts) : AmrCell() {
    build2D(verts);
}

AmrCell::AmrCell(const VerticesList& verts) : AmrCell() {
    build2D(verts);
}

AmrCell::AmrCell(const SqQuad& verts) : AmrCell() {
    build2D(verts);
}

AmrCell::AmrCell(const Cube& verts) : AmrCell() {
    build3D(verts);
}

void AmrCell::copy_to(AmrCell& cell) const {
    cell = *this;
}

double AmrCell::volume() const {
    return size * (dim < 3 ? size : size * size);
}

void AmrCell::set_undefined() {
    rank = -1;
}

bool AmrCell::is_actual() const {
    return rank >= 0;
}

bool AmrCell::is_undefined() const {
    return rank < 0;
}

/*
AmrCell::AmrCell(const std::array<Storage::Item, 4>& children) : dim(2) {
    SqQuad verts = {
            children[0].geom().vertices[iww(0, 0)],
            children[0].geom().vertices[iww(2, 0)],
            children[1].geom().vertices[iww(2, 0)],
            children[0].geom().vertices[iww(0, 2)],
            children[0].geom().vertices[iww(2, 2)],
            children[1].geom().vertices[iww(2, 2)],
            children[2].geom().vertices[iww(0, 2)],
            children[2].geom().vertices[iww(2, 2)],
            children[3].geom().vertices[iww(2, 2)]
    };
    build2D(verts);
}

AmrCell::AmrCell(const std::array<Storage::Item, 8>& children) : dim(3) {
    Cube verts = {
            children[iss(0, 0, 0)].geom().vertices[isw(0, 0, 0)],
            children[iss(1, 0, 0)].geom().vertices[isw(1, 0, 0)],
            children[iss(0, 1, 0)].geom().vertices[isw(0, 1, 0)],
            children[iss(1, 1, 0)].geom().vertices[isw(1, 1, 0)],
            children[iss(0, 0, 1)].geom().vertices[isw(0, 0, 1)],
            children[iss(1, 0, 1)].geom().vertices[isw(1, 0, 1)],
            children[iss(0, 1, 1)].geom().vertices[isw(0, 1, 1)],
            children[iss(1, 1, 1)].geom().vertices[isw(1, 1, 1)]
    };
    build3D(verts);
}
*/

} // namespace zephyr::geom