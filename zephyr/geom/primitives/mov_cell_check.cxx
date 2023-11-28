#include <fstream>
#include <iomanip>

#include <zephyr/geom/primitives/mov_cell.h>

namespace zephyr::geom {
/*
void PolyCell::print_info() const {
    std::cout << "\t\tcenter: " << coords.transpose() << "\n";
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
    for (int i = 0; i < 27; ++i) {
        if (vertices[i].is_actual()) {
            std::cout << "\t\t\t" << i << ": " << vertices[i].transpose() << "\n";
        }
    }

    std::cout << "\t\tcell.faces:\n";
    for (int i = 0; i < BFaces::max_size; ++i) {
        auto &face = faces[i];
        if (face.is_undefined()) continue;

        std::cout << "\t\t\t" << side_to_string(Side(i % 6)) << " face (" << i / 6 << "):\n";
        std::cout << "\t\t\t\tvertices:";
        for (int j = 0; j < (dim < 3 ? 2 : 4); ++j) {
            std::cout << " " << face.vertices[j];
        }
        std::cout << "\n";
        std::cout << "\t\t\t\tflag:       " << face.boundary << "\n";
        std::cout << "\t\t\t\tarea:       " << face.area << "\n";
        std::cout << "\t\t\t\tnormal:     " << face.normal.transpose() << "\n";
        std::cout << "\t\t\t\tadj.rank:   " << face.adjacent.rank << "\n";
        std::cout << "\t\t\t\tadj.index:  " << face.adjacent.index << "\n";
        std::cout << "\t\t\t\tadj.ghost:  " << face.adjacent.ghost << "\n";
    }
}

void PolyCell::visualize() const {
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
    SqQuad map;
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

    SqLine L(SqLine({vertices[0], vertices[3], vertices[6]}));
    SqLine R(SqLine({vertices[2], vertices[5], vertices[8]}));
    SqLine B(SqLine({vertices[0], vertices[1], vertices[2]}));
    SqLine T(SqLine({vertices[6], vertices[7], vertices[8]}));

    for (auto& sf: {L, R, B, T}) {
        file << "curve_Lx = spline(" << sf.v1.x() << ", " << sf.vc.x() << ", " << sf.v2.x() << ", unit)\n";
        file << "curve_Ly = spline(" << sf.v1.y() << ", " << sf.vc.y() << ", " << sf.v2.y() << ", unit)\n";
        file << "ax.plot(curve_Lx, curve_Ly, linestyle='dotted', color='green', linewidth=0.5)\n\n";
    }

    for (int i = 0; i < BFaces::max_size; ++i) {
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

int PolyCell::check_geometry() const {
    for (auto &face: faces) {
        if (face.is_undefined()) continue;

        Vector3d fc(0.0, 0.0, 0.0);
        for (int iv = 0; iv < VpF(dim); ++iv) {
            fc += vertices[face.vertices[iv]];
        }
        fc /= VpF(dim);

        // Нормаль внешняя
        if (face.normal.dot(fc - coords) < 0.0) {
            std::cout << "\tWrong normal direction (inside cell)\n";
            print_info();
            return -1;
        }

        // Вершины грани перечислены в правильном порядке
        if (dim > 2) {
            Vector3d v0 = vertices[face.vertices[0]];
            Vector3d v1 = vertices[face.vertices[1]];
            Vector3d v2 = vertices[face.vertices[2]];
            Vector3d v3 = vertices[face.vertices[3]];

            Vector3d n1 = (v2 - v1).cross(v0 - v1);
            Vector3d n2 = (v1 - v2).cross(v3 - v2);
            if (n1.dot(n2) < 0.0) {
                std::cout << "Wrong order of vertices on face\n";
                print_info();
                return -1;
            }
        }
    }
    return 0;
}

int PolyCell::check_base_face_orientation() const {
    if (dim == 2) {
        Vector3d nx1 = faces[Side::L].normal;
        Vector3d nx2 = faces[Side::R].normal;
        Vector3d ny1 = faces[Side::B].normal;
        Vector3d ny2 = faces[Side::T].normal;

        if (nx1.dot(nx2) > -0.8) {
            std::cout << "\tOpposite outward normals (left-right) are co-directed\n";
            print_info();
            return -1;
        }
        if (ny1.dot(ny2) > -0.8) {
            std::cout << "\tOpposite outward normals (bottom-top) are co-directed\n";
            print_info();
            return -1;
        }
        if (nx1.cross(ny1).z() < 0.0) {
            std::cout << "\tWrong face orientation (left-bottom)\n";
            print_info();
            return -1;
        }
        if (nx2.cross(ny2).z() < 0.0) {
            std::cout << "\tWrong face orientation (right-top)\n";
            print_info();
            return -1;
        }
    } else {
        Vector3d nx1 = faces[Side::L].normal;
        Vector3d nx2 = faces[Side::R].normal;
        Vector3d ny1 = faces[Side::B].normal;
        Vector3d ny2 = faces[Side::T].normal;
        Vector3d nz1 = faces[Side::X].normal;
        Vector3d nz2 = faces[Side::F].normal;

        if (nx1.dot(nx2) > -0.8) {
            std::cout << "\tOpposite outward normals (left-right) are co-directed\n";
            print_info();
            return -1;
        }
        if (ny1.dot(ny2) > -0.8) {
            std::cout << "\tOpposite outward normals (bottom-top) are co-directed\n";
            print_info();
            return -1;
        }
        if (nz1.dot(nz2) > -0.8) {
            std::cout << "\tOpposite outward normals (back-front) are co-directed\n";
            print_info();
            return -1;
        }
        if (nx1.dot(ny1.cross(nz1)) > 0.0) {
            std::cout << "\tWrong face orientation (left-bottom-back)\n";
            print_info();
            return -1;
        }
        if (nx2.dot(ny2.cross(nz2)) < 0.0) {
            std::cout << "\tWrong face orientation (right-top-front)\n";
            print_info();
            return -1;
        }
    }
    return 0;
}

int PolyCell::check_base_vertices_order() const {
    if (dim == 2) {
        SqQuad verts;
        for (int i = 0; i < 9; ++i) {
            verts[i] = vertices[i];
        }

        bool bad = false;

        // Индекс пересечения двух граней
        auto cross_face = [](const BFaces& faces, Side side1, Side side2) -> int {
            auto& face1 = faces[side1];
            auto& face2 = faces[side2];
            if (face1.is_undefined()) {
                return 100;
            }
            if (face2.is_undefined()) {
                return 100;
            }
            for (int i: {0, 1}) {
                for (int j: {0, 1}) {
                    if (face1.vertices[i] == face2.vertices[j]) {
                        return face1.vertices[i];
                    }
                }
            }
            return 100;
        };

        for (int i: {0, 1}) {
            for (int j: {0, 1}) {
                Vector3d a = vertices[iww(i + 1, j)] - vertices[iww(i, j)];
                Vector3d b = vertices[iww(i, j + 1)] - vertices[iww(i, j)];

                Vector3d c = vertices[iww(i, j + 1)] - vertices[iww(i + 1, j + 1)];
                Vector3d d = vertices[iww(i + 1, j)] - vertices[iww(i + 1, j + 1)];

                if (a.cross(b).z() < 0.0 ||
                    c.cross(d).z() < 0.0) {
                    bad = true;
                    break;
                }
            }
            if (bad) {
                break;
            }
        }

        // Пересечения граней по нужным вершинам
        if (cross_face(faces, Side::LEFT0, Side::BOTTOM0) != iww(0, 0)) {
            bad = true;
        }
        if (cross_face(faces, Side::LEFT0, Side::TOP) != iww(0, 2) &&
            cross_face(faces, Side::LEFT1, Side::TOP) != iww(0, 2)) {
            bad = true;
        }
        if (cross_face(faces, Side::RIGHT, Side::BOTTOM0) != iww(2, 0) &&
            cross_face(faces, Side::RIGHT, Side::BOTTOM1) != iww(2, 0)) {
            bad = true;
        }
        if (cross_face(faces, Side::RIGHT0, Side::TOP0) != iww(2, 2) &&
            cross_face(faces, Side::RIGHT0, Side::TOP1) != iww(2, 2) &&
            cross_face(faces, Side::RIGHT1, Side::TOP0) != iww(2, 2) &&
            cross_face(faces, Side::RIGHT1, Side::TOP1) != iww(2, 2)) {
            bad = true;
        }

        if (bad) {
            std::cout << "\tBad arrangement of vertices in cell\n";
            print_info();
            return -1;
        }
    } else {
        SqCube verts;
        for (int i = 0; i < 27; ++i) {
            verts[i] = vertices[i];
        }

        bool bad = false;

        // Индекс пересечения трех граней
        auto cross_face = [](const BFaces& faces, Side side1, Side side2, Side side3) -> int {
            auto& face1 = faces[side1];
            auto& face2 = faces[side2];
            auto& face3 = faces[side3];

            if (face1.is_undefined()) {
                return 100;
            }
            if (face2.is_undefined()) {
                return 100;
            }
            if (face3.is_undefined()) {
                return 100;
            }
            for (int i: {0, 1, 2, 3}) {
                for (int j: {0, 1, 2, 3}) {
                    for (int k: {0, 1, 2, 3}) {
                        if (face1.vertices[i] == face2.vertices[j] &&
                            face2.vertices[j] == face3.vertices[k]) {
                            return face1.vertices[i];
                        }
                    }
                }
            }
            return 100;
        };

        for (int i: {0, 1}) {
            for (int j: {0, 1}) {
                for (int k: {0, 1}) {
                    Vector3d a = vertices[iww(i + 1, j, k)] - vertices[iww(i, j, k)];
                    Vector3d b = vertices[iww(i, j + 1, k)] - vertices[iww(i, j, k)];
                    Vector3d c = vertices[iww(i, j, k + 1)] - vertices[iww(i, j, k)];

                    Vector3d A = vertices[iww(i + 1, j + 1, k + 1)] - vertices[iww(i, j + 1, k + 1)];
                    Vector3d B = vertices[iww(i + 1, j + 1, k + 1)] - vertices[iww(i + 1, j, k + 1)];
                    Vector3d C = vertices[iww(i + 1, j + 1, k + 1)] - vertices[iww(i + 1, j + 1, k)];

                    if (a.dot(b.cross(c)) < 0.0 || A.dot(B.cross(C)) < 0.0) {
                        bad = true;
                        break;
                    }
                }
                if (bad) {
                    break;
                }
            }
            if (bad) {
                break;
            }
        }

        // Пересечения граней по нужным вершинам ???

        if (bad) {
            std::cout << "\tBad arrangement of vertices in cell\n";
            print_info();
            return -1;
        }
    }

    return 0;
}

int PolyCell::check_complex_faces() const {
    if (dim == 2) {
        for (int s = 0; s < 4; ++s) {
            auto f1 = faces[s];
            auto f2 = faces[s + 6];
            if (f2.is_undefined()) {
                continue;
            }

            if (std::abs(f1.normal.dot(f2.normal) - 1.0) > 1.0e-2) {
                std::cout << "\tSubfaces are not co-directed (" + side_to_string(s) + " side)\n";
                print_info();
                return -1;
            }
            int count = 0;
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    if (f1.vertices[i] == f2.vertices[j]) {
                        ++count;
                    }
                }
            }
            if (count != 1) {
                std::cout << "\tStrange complex face\n";
                print_info();
                return -1;
            }
        }
    } else {
        for (int s = 0; s < 6; ++s) {
            auto f1 = faces[s];
            auto f2 = faces[s + 6];
            if (f2.is_undefined()) {
                continue;
            }

            auto f3 = faces[s + 12];
            auto f4 = faces[s + 18];
            if (f3.is_undefined() || f4.is_undefined()) {
                std::cout << "\tComplex 3D face (" + side_to_string(s) + " side) has less than 4 subfaces\n";
                print_info();
                return -1;
            }

            double d1 = std::abs(f1.normal.dot(f2.normal) - 1.0);
            double d2 = std::abs(f1.normal.dot(f3.normal) - 1.0);
            double d3 = std::abs(f1.normal.dot(f4.normal) - 1.0);
            if (d1 + d2 + d3 > 1.0e-5) {
                std::cout << "\tSubfaces are not co-directed (" + side_to_string(s) + " side)\n";
                print_info();
                return -1;
            }
        }
    }

    return 0;
}
*/
} // namespace zephyr::geom