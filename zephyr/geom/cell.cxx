#include <fstream>
#include <iomanip>

#include <zephyr/geom/maps.h>
#include <zephyr/geom/cell.h>

namespace zephyr { namespace geom {

void Cell::setup_vertices(const ShortList2D& vlist) {
    using topology::iv;
    using topology::iss;
    using topology::iws;
    using topology::iww;

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
    using topology::iv;
    using topology::isw;
    using topology::iws;
    using topology::iww;

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
    using topology::iv;
    using topology::iss;
    using topology::iws;
    using topology::iww;

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
    using topology::face_indices;

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
    using topology::face_indices;

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

    short undef = topology::undef_index();

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
    using topology::face_indices;

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
    using topology::iww;

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
    using topology::iss;
    using topology::isw;

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

void Cell::print_info() const {
    std::cout << "\t\tcenter: " << coords << "\n";
    std::cout << "\t\trank:   " << rank << "\n";
    std::cout << "\t\tindex:  " << index << "\n";
    std::cout << "\t\tsize:   " << size << "\n";
    std::cout << "\t\tvolume: " << volume() << "\n";
    std::cout << "\t\tflag:   " << flag << "\n";
    std::cout << "\t\tnext:   " << next << "\n";
    std::cout << "\t\tb_idx:  " << b_idx << "\n";
    std::cout << "\t\tlevel:  " << level << "\n";
    std::cout << "\t\tz_idx:  " << z_idx << "\n";

    std::cout << "\t\tcell.vertices:\n";
    for (int i = 0; i < Vertices::max_size; ++i) {
        if (vertices[i].is_actual()) {
            std::cout << "\t\t\t" << i << ": " << vertices[i] << "\n";
        }
    }

    std::cout << "\t\tcell.faces:\n";
    for (int i = 0; i < Faces::max_size; ++i) {
        auto &face = faces[i];
        if (face.is_undefined()) continue;

        std::cout << "\t\t\t" << side_to_string(Side(i % 6)) << " face (" << i / 6 << "):\n";
        std::cout << "\t\t\t\tvertices:";
        for (int j = 0; j < (dim < 3 ? 2 : 4); ++j) {
            std::cout << " " << face.vertices[j];
        }
        std::cout << "\n";
        std::cout << "\t\t\t\tflag:       " << boundary_to_string(face.boundary) << "\n";
        std::cout << "\t\t\t\tarea:       " << face.area << "\n";
        std::cout << "\t\t\t\tnormal:     " << face.normal << "\n";
        std::cout << "\t\t\t\tadj.rank:   " << face.adjacent.rank << "\n";
        std::cout << "\t\t\t\tadj.index:  " << face.adjacent.index << "\n";
        std::cout << "\t\t\t\tadj.ghost:  " << face.adjacent.ghost << "\n";
    }
}

void Cell::visualize() const {
    std::ofstream file("cell.py");

    file << std::scientific << std::setprecision(6);

    file << "#!/bin/python3\n";
    file << "import numpy as np\n";
    file << "import matplotlib.pyplot as plt\n\n\n";

    file << "def spline(v1, vc, v2, t):\n";
    file << "    return vc + 0.5*(v2 - v1)*t + 0.5*(v2 - 2*vc + v1)*t**2\n\n\n";

    file << "def plot_arrow(ax, x, y, vx, vy, color, width=0.5):\n";
    file << "    phi = 0.1*np.pi\n";
    file << "    ax.plot([x - (+vx*np.cos(phi) + vy*np.sin(phi)), x, x - (vx*np.cos(phi) - vy*np.sin(phi))],\n"
         << "            [y - (-vx*np.sin(phi) + vy*np.cos(phi)), y, y - (vx*np.sin(phi) + vy*np.cos(phi))],\n"
         << "            color=color)\n\n";

    // Основные точки
    LargeList2D map;
    for (int i = 0; i < 9; ++i) {
        map[i] = vertices[i];
    }

    file << "fig = plt.figure(dpi=150, figsize=(8, 8))\n";
    file << "ax = fig.add_subplot()\n\n";
    file << "ax.set_aspect('equal')\n";

    file << "ax.plot([";
    for (int i = 0; i < 8; ++i) {
        file << map[i].x() << ", ";
    }
    file << map[8].x() << "],\n";
    file << "        [";
    for (int i = 0; i < 8; ++i) {
        file << map[i].y() << ", ";
    }
    file << map[8].y() << "],\n";
    file << "        linestyle='none', color='orange', marker='o')\n\n";

    file << "ax.plot([" << coords.x() << "], [" << coords.y() << "], color='black', marker='x')\n\n";

    for (int i = 0; i < 9; ++i) {
        file << "ax.text(" << vertices[i].x() << ", " << vertices[i].y() << ", " << i << ")\n";
    }

    file << "unit = np.linspace(-1.0, 1.0, 101)\n\n";

    Mapping1D L(LargeList1D({vertices[0], vertices[3], vertices[6]}));
    Mapping1D R(LargeList1D({vertices[2], vertices[5], vertices[8]}));
    Mapping1D B(LargeList1D({vertices[0], vertices[1], vertices[2]}));
    Mapping1D T(LargeList1D({vertices[6], vertices[7], vertices[8]}));

    for (auto& sf: {L, R, B, T}) {
        file << "curve_Lx = spline(" << sf.v1.x() << ", " << sf.vc.x() << ", " << sf.v2.x() << ", unit)\n";
        file << "curve_Ly = spline(" << sf.v1.y() << ", " << sf.vc.y() << ", " << sf.v2.y() << ", unit)\n";
        file << "ax.plot(curve_Lx, curve_Ly, linestyle='dotted', color='green', linewidth=0.5)\n\n";
    }

    for (int i = 0; i < Faces::max_size; ++i) {
        auto &face = faces[i];
        if (face.is_undefined()) {
            continue;
        }

        double area = face.area;
        Vector3d normal = face.normal;

        Vector3d v1 = vertices[face.vertices[0]];
        Vector3d v2 = vertices[face.vertices[1]];
        Vector3d vc = (v1 + v2) / 2.0;

        file << "# face " << i << "\n";

        file << "ax.plot([" << v1.x() << ", " << v2.x() << "], ["
             << v1.y() << ", " << v2.y() << "], color='green', marker='.')\n";

        double a = 0.03;
        file << "plot_arrow(ax, " << vc.x() << ", " << vc.y() << ", "
             << a * (v2.x() - v1.x()) << ", " << a * (v2.y() - v1.y()) << ", color='green')\n";

        double b = 0.25 * area;
        file << "ax.plot([" << vc.x() << ", " << vc.x() + b * normal.x() << "], ["
             << vc.y() << ", " << vc.y() + b * normal.y() << "], color='red')\n";

        double c = 0.2 * b;
        file << "plot_arrow(ax, " << vc.x() + b * normal.x() << ", " << vc.y() + b * normal.y() << ", "
             << c * normal.x() << ", " << c * normal.y() << ", color='red')\n\n";
    }

    file << "\n";
    file << "fig.tight_layout()\n";
    file << "plt.show()\n";
}

} // namespace geom
} // namespace zephyr