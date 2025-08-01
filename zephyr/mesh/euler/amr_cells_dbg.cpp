#include <iostream>
#include <iomanip>
#include <fstream>

#include <zephyr/utils/mpi.h>
#include <zephyr/geom/geom.h>
#include <zephyr/mesh/euler/amr_cells.h>

using zephyr::utils::mpi;
using namespace zephyr::geom;

namespace zephyr::mesh {

inline int VpF(int dim) { return dim == 2 ? 2 : 4; }

void AmrCells::print_info(index_t ic) const {
    std::cout << "\t\tarr.ic: " << ic << "\n";
    std::cout << "\t\tcenter: " << center[ic].transpose() << "\n";
    std::cout << "\t\trank:   " << rank[ic] << "\n";
    std::cout << "\t\tindex:  " << index[ic] << "\n";
    std::cout << "\t\tsize:   " << linear_size(ic) << "\n";
    std::cout << "\t\tvolume: " << volume[ic] << "\n";
    std::cout << "\t\tflag:   " << flag[ic] << "\n";
    std::cout << "\t\tnext:   " << next[ic] << "\n";
    std::cout << "\t\tb_idx:  " << b_idx[ic] << "\n";
    std::cout << "\t\tlevel:  " << level[ic] << "\n";
    std::cout << "\t\tz_idx:  " << z_idx[ic] << "\n";

    std::cout << "\t\tcell.vertices:\n";
    for (index_t i: nodes_range(ic)) {
        if (!verts[i].hasNaN()) {
            std::cout << "\t\t\t" << i - node_begin[ic] << ": " << verts[i].transpose() << "\n";
        }
    }

    std::cout << "\t\tcell.faces:\n";
    for (index_t iface: faces_range(ic)) {
        if (faces.is_undefined(iface)) continue;

        std::cout << "\t\t\t" << side_to_string(iface - face_begin[ic], m_dim) << ":\n";
        std::cout << "\t\t\t\tvertices:";
        for (int j = 0; j < (m_dim < 3 ? 2 : 4); ++j) {
            std::cout << " " << faces.vertices[iface][j];
        }
        std::cout << "\n";
        std::cout << "\t\t\t\tflag:       " << faces.boundary[iface] << "\n";
        std::cout << "\t\t\t\tarea:       " << faces.area[iface] << "\n";
        std::cout << "\t\t\t\tnormal:     " << faces.normal[iface].transpose() << "\n";
        std::cout << "\t\t\t\tcenter:     " << faces.center[iface].transpose() << "\n";
        std::cout << "\t\t\t\tadj.rank:   " << faces.adjacent.rank[iface] << "\n";
        std::cout << "\t\t\t\tadj.index:  " << faces.adjacent.index[iface] << "\n";
        std::cout << "\t\t\t\tadj.alien:  " << faces.adjacent.alien[iface] << "\n";
        std::cout << "\t\t\t\tadj.basic:  " << faces.adjacent.basic[iface] << "\n";
    }
}

void AmrCells::visualize(index_t ic, std::string filename) const {
    if (filename.find(".py") == std::string::npos) {
        filename += ".py";
    }

    std::ofstream file(filename);

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
    SqQuad map = mapping<2>(ic);

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

    file << "ax.plot([" << center[ic].x() << "], [" << center[ic].y() << "], color='black', marker='x')\n\n";

    if (m_dim != 2) {
        throw std::runtime_error("Can't visualize 3D cell, sorry");
    }

    SqQuad vertices = mapping<2>(ic);

    for (int i = 0; i < 9; ++i) {
        file << "ax.text(" << vertices[i].x() << ", " << vertices[i].y() << ", " << i << ")\n";
    }

    file << "unit = np.linspace(-1.0, 1.0, 101)\n\n";

    SqLine L(SqLine({vertices[0], vertices[3], vertices[6]}));
    SqLine R(SqLine({vertices[2], vertices[5], vertices[8]}));
    SqLine B(SqLine({vertices[0], vertices[1], vertices[2]}));
    SqLine T(SqLine({vertices[6], vertices[7], vertices[8]}));

    for (auto& sf: {L, R, B, T}) {
        file << "curve_Lx = spline(" << sf(-1).x() << ", " << sf(0).x() << ", " << sf(+1).x() << ", unit)\n";
        file << "curve_Ly = spline(" << sf(-1).y() << ", " << sf(0).y() << ", " << sf(+1).y() << ", unit)\n";
        file << "ax.plot(curve_Lx, curve_Ly, linestyle='dotted', color='green', linewidth=0.5)\n\n";
    }

    for (index_t iface: faces_range(ic)) {
        if (faces.is_undefined(iface)) {
            continue;
        }

        double area = faces.area[iface];
        Vector3d normal = faces.normal[iface];

        Vector3d v1 = vertices[faces.vertices[iface][0]];
        Vector3d v2 = vertices[faces.vertices[iface][1]];
        Vector3d vc = faces.center[iface];

        file << "# face " << iface << "\n";

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

int AmrCells::check_geometry(index_t ic) const {
    for (index_t iface: faces_range(ic)) {
        if (faces.is_undefined(iface)) continue;

        Vector3d fc(0.0, 0.0, 0.0);
        for (int iv = 0; iv < VpF(m_dim); ++iv) {
            fc += verts[node_begin[ic] + faces.vertices[iface][iv]];
        }
        fc /= VpF(m_dim);

        // Нормаль внешняя
        if (faces.normal[iface].dot(fc - center[ic]) < 0.0) {
            std::cout << "\tWrong normal direction (inside cell)\n";
            print_info(ic);
            return -1;
        }

        // Вершины грани перечислены в правильном порядке
        if (m_dim > 2) {
            Vector3d v0 = verts[node_begin[ic] + faces.vertices[iface][0]];
            Vector3d v1 = verts[node_begin[ic] + faces.vertices[iface][1]];
            Vector3d v2 = verts[node_begin[ic] + faces.vertices[iface][2]];
            Vector3d v3 = verts[node_begin[ic] + faces.vertices[iface][3]];

            Vector3d n1 = (v2 - v1).cross(v0 - v1);
            Vector3d n2 = (v1 - v2).cross(v3 - v2);
            if (n1.dot(n2) < 0.0) {
                std::cout << "Wrong order of vertices on face\n";
                print_info(ic);
                return -1;
            }
        }
    }
    return 0;
}

int AmrCells::check_base_face_orientation(index_t ic) const {
    if (m_dim == 2) {
        Vector3d nx1 = faces.normal[face_begin[ic] + Side3D::L];
        Vector3d nx2 = faces.normal[face_begin[ic] + Side3D::R];
        Vector3d ny1 = faces.normal[face_begin[ic] + Side3D::B];
        Vector3d ny2 = faces.normal[face_begin[ic] + Side3D::T];

        if (nx1.dot(nx2) > -0.8) {
            std::cout << "\tOpposite outward normals (left-right) are co-directed\n";
            print_info(ic);
            return -1;
        }
        if (ny1.dot(ny2) > -0.8) {
            std::cout << "\tOpposite outward normals (bottom-top) are co-directed\n";
            print_info(ic);
            return -1;
        }
        if (nx1.cross(ny1).z() < 0.0) {
            std::cout << "\tWrong face orientation (left-bottom)\n";
            print_info(ic);
            return -1;
        }
        if (nx2.cross(ny2).z() < 0.0) {
            std::cout << "\tWrong face orientation (right-top)\n";
            print_info(ic);
            return -1;
        }
    } else {
        Vector3d nx1 = faces.normal[face_begin[ic] + Side3D::L];
        Vector3d nx2 = faces.normal[face_begin[ic] + Side3D::R];
        Vector3d ny1 = faces.normal[face_begin[ic] + Side3D::B];
        Vector3d ny2 = faces.normal[face_begin[ic] + Side3D::T];
        Vector3d nz1 = faces.normal[face_begin[ic] + Side3D::X];
        Vector3d nz2 = faces.normal[face_begin[ic] + Side3D::F];

        if (nx1.dot(nx2) > -0.8) {
            std::cout << "\tOpposite outward normals (left-right) are co-directed\n";
            print_info(ic);
            return -1;
        }
        if (ny1.dot(ny2) > -0.8) {
            std::cout << "\tOpposite outward normals (bottom-top) are co-directed\n";
            print_info(ic);
            return -1;
        }
        if (nz1.dot(nz2) > -0.8) {
            std::cout << "\tOpposite outward normals (back-front) are co-directed\n";
            print_info(ic);
            return -1;
        }
        if (nx1.dot(ny1.cross(nz1)) > 0.0) {
            std::cout << "\tWrong face orientation (left-bottom-back)\n";
            print_info(ic);
            return -1;
        }
        if (nx2.dot(ny2.cross(nz2)) < 0.0) {
            std::cout << "\tWrong face orientation (right-top-front)\n";
            print_info(ic);
            return -1;
        }
    }
    return 0;
}

int AmrCells::check_base_vertices_order(index_t ic) const {
    const double h = linear_size(ic);
    auto close = [h](Vector3d& x, Vector3d& y) -> bool {
        return (x - y).norm() < 1.0e-6 * h;
    };

    if (m_dim == 2) {
        SqQuad quad = mapping<2>(ic);

        bool bad = false;

        // Индекс пересечения двух граней
        auto cross_face = [this](index_t iface_base, Side2D side1, Side2D side2) -> int {
            auto iface = iface_base + side1;
            auto jface = iface_base + side2;
            if (faces.is_undefined(iface)) {
                return 100;
            }
            if (faces.is_undefined(jface)) {
                return 100;
            }
            for (int i: {0, 1}) {
                for (int j: {0, 1}) {
                    if (faces.vertices[iface][i] == faces.vertices[jface][j]) {
                        return faces.vertices[iface][i];
                    }
                }
            }
            return 100;
        };

        for (int i: {-1, 0}) {
            for (int j: {-1, 0}) {
                Vector3d a = quad(i + 1, j) - quad(i, j);
                Vector3d b = quad(i, j + 1) - quad(i, j);

                Vector3d c = quad(i, j + 1) - quad(i + 1, j + 1);
                Vector3d d = quad(i + 1, j) - quad(i + 1, j + 1);

                if (a.cross(b).z() < 0.0 || c.cross(d).z() < 0.0) {
                    bad = true;
                    break;
                }
            }
            if (bad) {
                break;
            }
        }

        // Пересечения граней по нужным вершинам
        if (cross_face(face_begin[ic], Side2D::LEFT, Side2D::BOTTOM), SqQuad::iss<-1, -1>()) {
            bad = true;
        }
        if (cross_face(face_begin[ic], Side2D::LEFT[0], Side2D::TOP) != SqQuad::iss<-1, +1>() &&
            cross_face(face_begin[ic], Side2D::LEFT[1], Side2D::TOP) != SqQuad::iss<-1, +1>()) {
            bad = true;
        }
        if (cross_face(face_begin[ic], Side2D::RIGHT, Side2D::BOTTOM[0]) != SqQuad::iss<+1, -1>() &&
            cross_face(face_begin[ic], Side2D::RIGHT, Side2D::BOTTOM[1]) != SqQuad::iss<+1, -1>()) {
            bad = true;
        }
        if (cross_face(face_begin[ic], Side2D::RIGHT[0], Side2D::TOP[0]) != SqQuad::iss<+1, +1>() &&
            cross_face(face_begin[ic], Side2D::RIGHT[0], Side2D::TOP[1]) != SqQuad::iss<+1, +1>() &&
            cross_face(face_begin[ic], Side2D::RIGHT[1], Side2D::TOP[0]) != SqQuad::iss<+1, +1>() &&
            cross_face(face_begin[ic], Side2D::RIGHT[1], Side2D::TOP[1]) != SqQuad::iss<+1, +1>()) {
            bad = true;
        }

        if (bad) {
            std::cout << "\tBad arrangement of vertices in cell (2D)\n";
            print_info(ic);
            return -1;
        }
    } else {
        SqCube cube = mapping<3>(ic);

        bool bad = false;

        throw std::runtime_error("Can't check 3D cell");

        /*
        // Индекс пересечения трех граней
        auto cross_face = [](const BFaces& faces, Side3D side1, Side3D side2, Side3D side3) -> int {
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
         */

        for (int i: {-1, 0}) {
            for (int j: {-1, 0}) {
                for (int k: {-1, 0}) {
                    Vector3d a = cube(i + 1, j, k) - cube(i, j, k);
                    Vector3d b = cube(i, j + 1, k) - cube(i, j, k);
                    Vector3d c = cube(i, j, k + 1) - cube(i, j, k);

                    Vector3d A = cube(i + 1, j + 1, k + 1) - cube(i, j + 1, k + 1);
                    Vector3d B = cube(i + 1, j + 1, k + 1) - cube(i + 1, j, k + 1);
                    Vector3d C = cube(i + 1, j + 1, k + 1) - cube(i + 1, j + 1, k);

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
            std::cout << "\tBad arrangement of vertices in cell (3D)\n";
            print_info(ic);
            return -1;
        }
    }

    return 0;
}

int AmrCells::check_complex_faces(index_t ic) const {
    if (m_dim == 2) {
        for (Side2D side: Side2D::items()) {
            auto iface1 = face_begin[ic] + side;
            auto iface2 = face_begin[ic] + side[1];
            if (faces.is_undefined(iface2)) {
                continue;
            }

            if (std::abs(faces.normal[iface1].dot(faces.normal[iface2]) - 1.0) > 1.0e-2) {
                std::cout << "\tSubfaces are not co-directed (" << side_to_string(side, m_dim) << ")\n";
                print_info(ic);
                return -1;
            }
            int count = 0;
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    if (faces.vertices[iface1][i] == faces.vertices[iface2][j]) {
                        ++count;
                    }
                }
            }
            if (count != 1) {
                std::cout << "\tStrange complex face\n";
                print_info(ic);
                return -1;
            }
        }
    } else {
        for (Side3D side: Side3D::items()) {
            auto iface1 = face_begin[ic] + side;
            auto iface2 = face_begin[ic] + side[1];
            if (faces.is_undefined(iface2)) {
                continue;
            }

            auto iface3 = face_begin[ic] + side[2];
            auto iface4 = face_begin[ic] + side[3];
            if (faces.is_undefined(iface3) || faces.is_undefined(iface4)) {
                std::cout << "\tComplex 3D face (" + side_to_string(side, m_dim) + " side) has less than 4 subfaces\n";
                print_info(ic);
                return -1;
            }

            double d1 = std::abs(faces.normal[iface1].dot(faces.normal[iface2]) - 1.0);
            double d2 = std::abs(faces.normal[iface1].dot(faces.normal[iface3]) - 1.0);
            double d3 = std::abs(faces.normal[iface1].dot(faces.normal[iface4]) - 1.0);
            if (d1 + d2 + d3 > 1.0e-5) {
                std::cout << "\tSubfaces are not co-directed (" + side_to_string(side, m_dim) << ")\n";
                print_info(ic);
                return -1;
            }
        }
    }

    return 0;
}

int AmrCells::check_connectivity(index_t ic) const {
    AmrCells aliens;
    return check_connectivity(ic, aliens);
}

int AmrCells::check_connectivity(index_t ic, const AmrCells& aliens) const {
    if (ic >= m_size) {
        throw std::runtime_error("Данная проверка только для локальных ячеек!");
    }

    if (is_undefined(ic)) return 0;

    // Через обычные грани существуют соседи
    for (index_t iface: faces_range(ic)) {
        if (faces.is_undefined(iface)) {
            // Проверим обнуление параметров
            if (faces.adjacent.rank[iface] >= 0) {
                std::cout << "\tUndefined face, adjacent.rank >= 0\n";
                print_info(ic);
                return -1;
            }
            if (faces.adjacent.index[iface] >= 0) {
                std::cout << "\tUndefined face, adjacent.index >= 0\n";
                print_info(ic);
                return -1;
            }
            if (faces.adjacent.alien[iface] >= 0) {
                std::cout << "\tUndefined face, adjacent.alien >= 0\n";
                print_info(ic);
                return -1;
            }
            continue;
        }

        // Массив смежности
        auto& adj = faces.adjacent;

        std::string f_name = face_name(ic, iface);

        // Должна указывать на ячейку
        if (faces.adjacent.basic[iface] != ic) {
            std::cout << "\tFace " << f_name << " should point to origin cell\n";
            print_info(ic);
            return -1;
        }

        // Простая граничная грань
        if (faces.is_boundary(iface)) {
            // Грань должна ссылаться на саму ячейку
            if (adj.rank[iface] != rank[ic] || adj.index[iface] != ic || adj.alien[iface] >= 0) {
                std::cout << "\tBoundary face should point to origin cell\n";
                print_info(ic);
                return -1;
            }
            continue;
        }

        // Отрицательный ранг, почему?
        if (adj.rank[iface] < 0) {
            std::cout << "\tadjacent.rank < 0\n";
            print_info(ic);
            return -1;
        }

        if (adj.index[iface] < 0) {
            std::cout << "\tadjacent.index out of range #1\n";
            print_info(ic);
            return -1;
        }

        if (adj.rank[iface] == mpi::rank()) {
            // Локальный сосед
            if (adj.index[iface] >= m_size) {
                std::cout << "\tadjacent.index out of range #1\n";
                print_info(ic);
                return -1;
            }
            if (adj.alien[iface] >= 0) {
                std::cout << "\tadjacent.alien >= 0 for local cell\n";
                print_info(ic);
                return -1;
            }
            if (adj.index[iface] == ic) {
                std::cout << "\tSelf reference\n";
                print_info(ic);
                return -1;
            }
        }
        else {
            // Удаленная ячейка
            if (adj.alien[iface] < 0 || adj.alien[iface] >= aliens.size()) {
                std::cout << "\t" + f_name + ": adjacent.alien out of range for remote cell\n";
                std::cout << "\t\taliens.size: " << aliens.size() << "\n";
                print_info(ic);
                return -1;
            }
        }

        // Сосед может быть из aliens
        auto [neibs, jc] = adj.get_neib(iface, *this, aliens);

        Vector3d fc = faces.center[iface];

        // Сосед должен иметь точно такую же грань, но с противоположной нормалью,
        // и ссылаться на текущую ячейку
        int counter = 0;

        for (index_t jface: neibs.faces_range(jc)) {
            if (neibs.faces.is_undefined(jface)) {
                continue;
            }

            // простая граничная грань
            if (neibs.faces.is_boundary(jface)) {
                continue;
            }

            Vector3d nfc = neibs.faces.center[jface];
            if ((fc - nfc).norm() > 1.0e-6 * linear_size(ic)) {
                continue;
            }

            // мы нашли соответствующую грань соседа
            ++counter;

            // нормали противоположны
            if (std::abs(faces.normal[iface].dot(neibs.faces.normal[jface]) + 1.0) > 1.0e-6) {
                std::cout << "\tOpposite faces have not opposite normals\n";
                std::cout << "\tCurrent cell:\n";
                print_info(ic);
                std::cout << "\tNeighbor:\n";
                neibs.print_info(jc);
                return -1;
            }
            // площади совпадают
            if (std::abs(faces.area[iface] - neibs.faces.area[jface]) > 1.0e-6 * linear_size(ic)) {
                std::cout << "\tOpposite faces have different area\n";
                std::cout << "\tCurrent cell:\n";
                print_info(ic);
                std::cout << "\tNeighbor:\n";
                neibs.print_info(jc);
                return -1;
            }
            // Указывает на исходную ячейку
            if (neibs.faces.adjacent.index[jface] != ic) {
                std::cout << "\tWrong connection (index != ic). " << f_name << " face\n";
                std::cout << "\tCurrent cell:\n";
                print_info(ic);
                std::cout << "\tNeighbor:\n";
                neibs.print_info(jc);
                return -1;
            }
            // Ранг смежной (исходная) больше нуля
            if (neibs.faces.adjacent.rank[jface] < 0) {
                std::cout << "\tWrong adjacent (rank < 0)\n";
                std::cout << "\tCurrent cell:\n";
                print_info(ic);
                std::cout << "\tNeighbor:\n";
                neibs.print_info(jc);
                return -1;
            }
            // Ранг смежной (исходная) равен рангу процесса
            if (neibs.faces.adjacent.rank[jface] != mpi::rank()) {
                std::cout << "\tWrong adjacent (rank)\n";
                std::cout << "\tCurrent cell:\n";
                print_info(ic);
                std::cout << "\tNeighbor:\n";
                neibs.print_info(jc);
                return -1;
            }
        }
        if (counter < 1) {
            std::cout << "\tHas no neighbor across ordinary " << f_name << "\n";
            print_info(ic);
            std::cout << "\tNeighbor:\n";
            neibs.print_info(jc);
            return -1;
        }
        if (counter > 1) {
            std::cout << "\tMore than one neighbor across ordinary face\n";
            print_info(ic);
            return -1;
        }
    }

    return 0;
}

} // namespace zephyr::mesh