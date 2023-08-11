#include <zephyr/geom/cell.h>

namespace zephyr { namespace geom {

void Cell::setup_vertices(const ShortList2D& vlist) {
    vertices.set_undefined();

    // Базовые вершины
    vertices[iv(0, 0)] = vlist[iss(0, 0)];
    vertices[iv(0, 1)] = vlist[iss(0, 1)];
    vertices[iv(1, 0)] = vlist[iss(1, 0)];
    vertices[iv(1, 1)] = vlist[iss(1, 1)];

    // Вершины на гранях
    vertices[iww(0, 1)] = (vlist[iws(0, 0)] + vlist[iws(0, 2)]) / 2.0;
    vertices[iww(2, 1)] = (vlist[iws(2, 0)] + vlist[iws(2, 2)]) / 2.0;
    vertices[iww(1, 0)] = (vlist[iws(0, 0)] + vlist[iws(2, 0)]) / 2.0;
    vertices[iww(1, 2)] = (vlist[iws(0, 2)] + vlist[iws(2, 2)]) / 2.0;

    // Вершина в центре
    vertices[iww(1, 1)] = (vlist[iws(0, 0)] + vlist[iws(2, 0)] +
                           vlist[iws(0, 2)] + vlist[iws(2, 2)]) / 4.0;
}

void Cell::setup_vertices(const LargeList2D& vlist) {
    vertices.set_undefined();

    // Базовые вершины
    vertices[iv(0, 0)] = vlist[isw(0, 0)];
    vertices[iv(0, 1)] = vlist[isw(0, 1)];
    vertices[iv(1, 0)] = vlist[isw(1, 0)];
    vertices[iv(1, 1)] = vlist[isw(1, 1)];

    // Вершины на гранях
    vertices[iww(0, 1)] = vlist[iww(0, 1)];
    vertices[iww(2, 1)] = vlist[iww(2, 1)];
    vertices[iww(1, 0)] = vlist[iww(1, 0)];
    vertices[iww(1, 2)] = vlist[iww(1, 2)];

    // Вершина в центре
    vertices[iww(1, 1)] = vlist[iww(1, 1)];
}

void Cell::setup_vertices(const ShortList3D& vlist) {
    vertices.set_undefined();

    // Скопировать базовые вершины
    vertices[iv(0, 0, 0)] = vlist[iss(0, 0, 0)];
    vertices[iv(1, 0, 0)] = vlist[iss(1, 0, 0)];
    vertices[iv(0, 1, 0)] = vlist[iss(0, 1, 0)];
    vertices[iv(1, 1, 0)] = vlist[iss(1, 1, 0)];
    vertices[iv(0, 0, 1)] = vlist[iss(0, 0, 1)];
    vertices[iv(1, 0, 1)] = vlist[iss(1, 0, 1)];
    vertices[iv(0, 1, 1)] = vlist[iss(0, 1, 1)];
    vertices[iv(1, 1, 1)] = vlist[iss(1, 1, 1)];

    // Вершины на ребрах
    vertices[iww(0, 1, 0)] = (vlist[iws(0, 0, 0)] + vlist[iws(0, 2, 0)]) / 2.0;
    vertices[iww(2, 1, 0)] = (vlist[iws(2, 0, 0)] + vlist[iws(2, 2, 0)]) / 2.0;
    vertices[iww(1, 0, 0)] = (vlist[iws(0, 0, 0)] + vlist[iws(2, 0, 0)]) / 2.0;
    vertices[iww(1, 2, 0)] = (vlist[iws(0, 2, 0)] + vlist[iws(2, 2, 0)]) / 2.0;
    vertices[iww(0, 0, 1)] = (vlist[iws(0, 0, 0)] + vlist[iws(0, 0, 2)]) / 2.0;
    vertices[iww(2, 0, 1)] = (vlist[iws(2, 0, 0)] + vlist[iws(2, 0, 2)]) / 2.0;
    vertices[iww(0, 2, 1)] = (vlist[iws(0, 2, 0)] + vlist[iws(0, 2, 2)]) / 2.0;
    vertices[iww(2, 2, 1)] = (vlist[iws(2, 2, 0)] + vlist[iws(2, 2, 2)]) / 2.0;
    vertices[iww(0, 1, 2)] = (vlist[iws(0, 0, 2)] + vlist[iws(0, 2, 2)]) / 2.0;
    vertices[iww(2, 1, 2)] = (vlist[iws(2, 0, 2)] + vlist[iws(2, 2, 2)]) / 2.0;
    vertices[iww(1, 0, 2)] = (vlist[iws(0, 0, 2)] + vlist[iws(2, 0, 2)]) / 2.0;
    vertices[iww(1, 2, 2)] = (vlist[iws(0, 2, 2)] + vlist[iws(2, 2, 2)]) / 2.0;

    // Вершны на гранях
    vertices[iww(0, 1, 1)] = (vlist[iws(0, 0, 0)] + vlist[iws(0, 2, 0)] +
                                   vlist[iws(0, 0, 2)] + vlist[iws(0, 2, 2)]) / 4.0;
    vertices[iww(2, 1, 1)] = (vlist[iws(2, 0, 0)] + vlist[iws(2, 2, 0)] +
                                   vlist[iws(2, 0, 2)] + vlist[iws(2, 2, 2)]) / 4.0;
    vertices[iww(1, 0, 1)] = (vlist[iws(0, 0, 0)] + vlist[iws(2, 0, 0)] +
                                   vlist[iws(0, 0, 2)] + vlist[iws(2, 0, 2)]) / 4.0;
    vertices[iww(1, 2, 1)] = (vlist[iws(0, 2, 0)] + vlist[iws(2, 2, 0)] +
                                   vlist[iws(0, 2, 2)] + vlist[iws(2, 2, 2)]) / 4.0;
    vertices[iww(1, 1, 0)] = (vlist[iws(0, 0, 0)] + vlist[iws(2, 0, 0)] +
                                   vlist[iws(0, 2, 0)] + vlist[iws(2, 2, 0)]) / 4.0;
    vertices[iww(1, 1, 2)] = (vlist[iws(0, 0, 2)] + vlist[iws(2, 0, 2)] +
                                   vlist[iws(0, 2, 2)] + vlist[iws(2, 2, 2)]) / 4.0;

    // Вершина в центре
    vertices[iww(1, 1, 1)] = (vlist[iws(0, 0, 0)] + vlist[iws(0, 2, 0)] +
                                   vlist[iws(2, 0, 0)] + vlist[iws(2, 2, 0)] +
                                   vlist[iws(0, 0, 2)] + vlist[iws(0, 2, 2)] +
                                   vlist[iws(2, 0, 2)] + vlist[iws(2, 2, 2)]) / 8.0;
}

void Cell::build2D(const ShortList2D& verts) {
    dim = 2;

    double area = geom::area(verts);
    Vector3d centroid = geom::centroid(verts, area);

    coords = centroid;
    size = std::sqrt(area);

    setup_vertices(verts);

    faces.set_undefined();

    faces[Side::L].vertices = face_indices<2, Side::L>();
    faces[Side::R].vertices = face_indices<2, Side::R>();
    faces[Side::B].vertices = face_indices<2, Side::B>();
    faces[Side::T].vertices = face_indices<2, Side::T>();

    for (int i = 0; i < 4; ++i) {
        ShortList1D vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]]
        };

        faces[i].area = geom::length(vs);
        faces[i].normal = geom::normal(vs, centroid);
        faces[i].boundary = FaceFlag::ORDINARY;
    }
}

void Cell::build2D(const LargeList2D& verts) {
    dim = 2;

    double area = geom::area(verts);
    Vector3d centroid = geom::centroid(verts, area);

    coords = centroid;
    size = std::sqrt(area);

    setup_vertices(verts);

    faces.set_undefined();

    faces[Side::L].vertices = face_indices<2, Side::L>();
    faces[Side::R].vertices = face_indices<2, Side::R>();
    faces[Side::B].vertices = face_indices<2, Side::B>();
    faces[Side::T].vertices = face_indices<2, Side::T>();

    for (int i = 0; i < 4; ++i) {
        ShortList1D vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]]
        };

        faces[i].area = geom::length(vs);
        faces[i].normal = geom::normal(vs, centroid);
        faces[i].boundary = FaceFlag::ORDINARY;
    }
}

void Cell::build2D(const VerticesList& verts) {
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

        ShortList1D vs = {verts[i], verts[j]};

        faces[i].area = geom::length(vs);
        faces[i].normal = geom::normal(vs, centroid);
        faces[i].boundary = FaceFlag::ORDINARY;
    }
}

void Cell::build3D(const ShortList3D& verts) {
    dim = 3;

    double volume = geom::volume(verts);
    size = std::cbrt(volume);

    Vector3d centroid = geom::centroid(verts, volume);
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
        ShortList2D vs = {
                vertices[faces[i].vertices[0]],
                vertices[faces[i].vertices[1]],
                vertices[faces[i].vertices[2]],
                vertices[faces[i].vertices[3]]
        };

        faces[i].area = geom::area(vs);
        faces[i].normal = geom::normal(vs, centroid);
        faces[i].boundary = FaceFlag::ORDINARY;
    }
}

Cell::Cell() :
    rank(0),
    index(0),
    b_idx(0),
    z_idx(0),
    next(0),
    level(0),
    flag(0),
    dim(0),
    size(0.0) {

}

Cell::Cell(const ShortList2D& verts) : Cell() {
    build2D(verts);
}

Cell::Cell(const VerticesList& verts) : Cell() {
    build2D(verts);
}

Cell::Cell(const LargeList2D& verts) : Cell() {
    build2D(verts);
}

Cell::Cell(const ShortList3D& verts) : Cell() {
    build3D(verts);
}

void Cell::copy_to(Cell& cell) const {
    cell = *this;
}

double Cell::volume() const {
    return size * (dim < 3 ? size : size * size);
}

void Cell::set_undefined() {
    rank = -1;
}

bool Cell::is_actual() const {
    return rank >= 0;
}

bool Cell::is_undefined() const {
    return rank < 0;
}

/*
Cell::Cell(const std::array<Storage::Item, 4>& children) : dim(2) {
    LargeList2D vs = {
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
    build2D(vs);
}

Cell::Cell(const std::array<Storage::Item, 8>& children) : dim(3) {
    ShortList3D vs = {
            children[iss(0, 0, 0)].geom().vertices[isw(0, 0, 0)],
            children[iss(1, 0, 0)].geom().vertices[isw(1, 0, 0)],
            children[iss(0, 1, 0)].geom().vertices[isw(0, 1, 0)],
            children[iss(1, 1, 0)].geom().vertices[isw(1, 1, 0)],
            children[iss(0, 0, 1)].geom().vertices[isw(0, 0, 1)],
            children[iss(1, 0, 1)].geom().vertices[isw(1, 0, 1)],
            children[iss(0, 1, 1)].geom().vertices[isw(0, 1, 1)],
            children[iss(1, 1, 1)].geom().vertices[isw(1, 1, 1)]
    };
    build3D(vs);
}
*/

} // namespace geom
} // namespace zephyr