/*
#include <fstream>

#include <zephyr/mesh/amr2/check.h>
#include <zephyr/math/geom/maps.h>
#include <zephyr/math/geom/cell.h>

using namespace ::zephyr::data;

using amrData;

namespace zephyr { namespace mesh { namespace amr2 {

void test_storage(AmrStorage &cells) {
    bool nice = cells.has_type(coords) && cells.has_type(vertices) &&
                cells.has_type(faces) && cells.has_type(size) &&
                cells.has_type(element) && cells.has_type(amrData);
    if (!nice) {
        throw std::runtime_error("AmrStorage should contain following types for refinement: "
                                 "coords, vertices, faces, size, element, amrData");
    }
    if (cells.empty()) {
        throw std::runtime_error("Empty AmrStorage for refiner");
    }
}

void print_cell_info(AmrStorage::Item cell) {
    using Vector3d;
    using coords;
    using vertices;
    using faces;
    using element;
    using side;

    auto dim = cell[element].dimension;

    std::cout << "\t\tcell.center: " << (Vector3d &) cell[coords] << "\n";
    std::cout << "\t\tuid: " << cell[uid].value << "\n";
    std::cout << "\t\tsize: " << cell[size].value << "\n";
    std::cout << "\t\tvolume: " << std::pow(cell[size].value, dim) << "\n";
    std::cout << "\t\tflag: " << cell.flag << "\n";
    std::cout << "\t\tnext: " << cell.next << "\n";
    std::cout << "\t\tbase_id: " << cell.b_idx << "\n";
    std::cout << "\t\tlevel: " << cell.level << "\n";
    std::cout << "\t\tz: " << cell.z() << "\n";
    std::cout << "\t\tcell.vertices:\n";
    for (int i = 0; i < cell[vertices].size(); ++i) {
        std::cout << "\t\t\t" << i << ": " << (Vector3d &) cell.vertices[i] << "\n";
    }
    std::cout << "\t\tcell.faces:\n";
    for (int i = 0; i < Faces::max_size; ++i) {
        auto &face = cell.faces[i];
        if (face.is_undefined()) continue;

        std::cout << "\t\t\t" << side_to_string(side(i % 6)) << " face (" << i / 6 << "):\n";
        std::cout << "\t\t\t\tvertices:";
        for (int j = 0; j < VpF(dim); ++j) {
            std::cout << " " << face.vertices[j];
        }
        std::cout << "\n";
        std::cout << "\t\t\t\tflag:   " << boundary_to_string(face.boundary) << "\n";
        std::cout << "\t\t\t\tarea:   " << face.area << "\n";
        std::cout << "\t\t\t\tnormal: " << (Vector3d &) face.normal << "\n";
        std::cout << "\t\t\t\tadj.rank: " << face.adjacent.rank << "\n";
        std::cout << "\t\t\t\tadj.index:  " << face.adjacent.index << "\n";
        std::cout << "\t\t\t\tadj.ghost:  " << face.adjacent.ghost << "\n";
    }
}

void print_cell_info(AmrStorage& locals, AmrStorage& aliens, size_t ic) {
    auto cell = locals[ic];
    print_cell_info(cell);
    std::cout << "\tAll neighbors of cell:\n";
    for (int i = 0; i < Faces::max_size; ++i) {
        auto &face = cell.faces[i];
        if (face.is_undefined() or face.is_boundary()) continue;

        std::cout << "\tNeighbor through the " << side_to_string(side(i % 6)) << " face (" << i / 6 << "):\n";

        if (face.adjacent.ghost > std::numeric_limits<int>::max()) {
            // Локальная ячейка
            if (face.adjacent.index >= locals.size()) {
                std::cout << "print_cell_info: Wrong connection #1 (It's acceptable for some intermediate refinement stages)";
            }
            else {
                auto neib = locals[face.adjacent.index];
                print_cell_info(neib);
            }
        }
        else {
            // Удаленная ячейка
            if (face.adjacent.ghost >= aliens.size()) {
                std::cout << "print_cell_info: Wrong connection #2 (It's acceptable for some intermediate refinement stages)";
            }
            else {
                auto neib = aliens[face.adjacent.ghost];
                print_cell_info(neib);
            }
        }
    }
}

void visualize_cell(AmrStorage::Item cell) {
    using zephyr::math::geom::SqLine;
    using zephyr::math::geom::Mapping2D;
    using zephyr::math::geom::LargeList1D;
    using zephyr::math::geom::LargeList2D;

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
    auto& faces = cell[faces];
    LargeList2D vertices;
    for (int i = 0; i < 9; ++i) {
        vertices[i] = (Vector3d &) cell.vertices[i];
    }

    file << "fig = plt.figure(dpi=150, figsize=(8, 8))\n";
    file << "ax = fig.add_subplot()\n\n";
    file << "ax.set_aspect('equal')\n";

    file << "ax.plot([";
    for (int i = 0; i < 8; ++i) {
        file << vertices[i].x << ", ";
    }
    file << vertices[8].x << "],\n";
    file << "        [";
    for (int i = 0; i < 8; ++i) {
        file << vertices[i].y << ", ";
    }
    file << vertices[8].y << "],\n";
    file << "        linestyle='none', color='orange', marker='o')\n\n";

    file << "ax.plot([" << cell[coords].x << "], [" << cell[coords].y << "], color='black', marker='x')\n\n";

    for (int i = 0; i < 9; ++i) {
        file << "ax.text(" << vertices[i].x << ", " << vertices[i].y << ", " << i << ")\n";
    }

    file << "unit = np.linspace(-1.0, 1.0, 101)\n\n";

    SqLine L(LargeList1D({vertices[0], vertices[3], vertices[6]}));
    SqLine R(LargeList1D({vertices[2], vertices[5], vertices[8]}));
    SqLine B(LargeList1D({vertices[0], vertices[1], vertices[2]}));
    SqLine T(LargeList1D({vertices[6], vertices[7], vertices[8]}));

    for (auto& map: {L, R, B, T}) {
        file << "curve_Lx = spline(" << map.v1.x << ", " << map.vc.x << ", " << map.v2.x << ", unit)\n";
        file << "curve_Ly = spline(" << map.v1.y << ", " << map.vc.y << ", " << map.v2.y << ", unit)\n";
        file << "ax.plot(curve_Lx, curve_Ly, linestyle='dotted', color='green', linewidth=0.5)\n\n";
    }

    for (int i = 0; i < Faces::max_size; ++i) {
        auto &face = faces[i];
        if (face.is_undefined()) {
            continue;
        }

        double area = face.area;
        Vector3d normal = (Vector3d &) face.normal;

        Vector3d v1 = vertices[face.vertices[0]];
        Vector3d v2 = vertices[face.vertices[1]];
        Vector3d vc = (v1 + v2) / 2.0;

        file << "# face " << i << "\n";

        file << "ax.plot([" << v1.x << ", " << v2.x << "], ["
             << v1.y << ", " << v2.y << "], color='green', marker='.')\n";

        double a = 0.03;
        file << "plot_arrow(ax, " << vc.x << ", " << vc.y << ", "
             << a * (v2.x - v1.x) << ", " << a * (v2.y - v1.y) << ", color='green')\n";

        double b = 0.25 * area;
        file << "ax.plot([" << vc.x << ", " << vc.x + b * normal.x << "], ["
             << vc.y << ", " << vc.y + b * normal.y << "], color='red')\n";

        double c = 0.2 * b;
        file << "plot_arrow(ax, " << vc.x + b * normal.x << ", " << vc.y + b * normal.y << ", "
             << c * normal.x << ", " << c * normal.y << ", color='red')\n\n";
    }

    Mapping2D map2D(vertices);

    file << "\n";
    file << "fig.tight_layout()\n";
    file << "plt.show()\n";
}

int check_geometry(zephyr::data::AmrStorage::Item cell) {
    using faces;
    using vertices;
    using Vector3d;
    using element;

    auto dim = cell[element].dimension;

    for (auto &face: cell.faces) {
        if (face.is_undefined()) continue;

        Vector3d fc(0.0, 0.0, 0.0);
        for (int iv = 0; iv < VpF(dim); ++iv) {
            fc += (Vector3d &) cell.vertices[face.vertices[iv]];
        }
        fc /= VpF(dim);

        // Нормаль внешняя
        if (scalar_product(face.normal, fc - (Vector3d &) cell[coords]) < 0.0) {
            std::cout << "\tWrong normal direction (inside cell)\n";
            print_cell_info(cell);
            return -1;
        }

        // Вершины грани перечислены в правильном порядке
        if (dim > 2) {
            Vector3d v0 = (Vector3d&)cell.vertices[face.vertices[0]];
            Vector3d v1 = (Vector3d&)cell.vertices[face.vertices[1]];
            Vector3d v2 = (Vector3d&)cell.vertices[face.vertices[2]];
            Vector3d v3 = (Vector3d&)cell.vertices[face.vertices[3]];

            Vector3d n1 = cross_product(v2 - v1, v0 - v1);
            Vector3d n2 = cross_product(v1 - v2, v3 - v2);
            if (scalar_product(n1, n2) < 0.0) {
                std::cout << "Wrong order of vertices on face\n";
                print_cell_info(cell);
                return -1;
            }
        }
    }
    return 0;
}

int check_base_face_orientation(AmrStorage::Item &cell) {
    using Vector3d;
    using faces;
    using vertices;
    using side;

    if (cell[element].dimension == 2) {
        Vector3d nx1 = (Vector3d &) cell.faces[Side::LEFT].normal;
        Vector3d nx2 = (Vector3d &) cell.faces[Side::RIGHT].normal;
        Vector3d ny1 = (Vector3d &) cell.faces[Side::BOTTOM].normal;
        Vector3d ny2 = (Vector3d &) cell.faces[Side::TOP].normal;

        if (scalar_product(nx1, nx2) > -0.8) {
            std::cout << "\tOpposite outward normals (left-right) are co-directed\n";
            print_cell_info(cell);
            return -1;
        }
        if (scalar_product(ny1, ny2) > -0.8) {
            std::cout << "\tOpposite outward normals (bottom-top) are co-directed\n";
            print_cell_info(cell);
            return -1;
        }
        if (nx1.x * ny1.y - nx1.y * ny1.x < 0.0) {
            std::cout << "\tWrong face orientation (left-bottom)\n";
            print_cell_info(cell);
            return -1;
        }
        if (nx2.x * ny2.y - nx2.y * ny2.x < 0.0) {
            std::cout << "\tWrong face orientation (right-top)\n";
            print_cell_info(cell);
            return -1;
        }
    } else {
        Vector3d nx1 = (Vector3d &) cell.faces[Side::LEFT].normal;
        Vector3d nx2 = (Vector3d &) cell.faces[Side::RIGHT].normal;
        Vector3d ny1 = (Vector3d &) cell.faces[Side::BOTTOM].normal;
        Vector3d ny2 = (Vector3d &) cell.faces[Side::TOP].normal;
        Vector3d nz1 = (Vector3d &) cell.faces[Side::BACK].normal;
        Vector3d nz2 = (Vector3d &) cell.faces[Side::FRONT].normal;

        if (scalar_product(nx1, nx2) > -0.8) {
            std::cout << "\tOpposite outward normals (left-right) are co-directed\n";
            print_cell_info(cell);
            return -1;
        }
        if (scalar_product(ny1, ny2) > -0.8) {
            std::cout << "\tOpposite outward normals (bottom-top) are co-directed\n";
            print_cell_info(cell);
            return -1;
        }
        if (scalar_product(nz1, nz2) > -0.8) {
            std::cout << "\tOpposite outward normals (back-front) are co-directed\n";
            print_cell_info(cell);
            return -1;
        }
        if (triple_product(nx1, ny1, nz1) > 0.0) {
            std::cout << "\tWrong face orientation (left-bottom-back)\n";
            print_cell_info(cell);
            return -1;
        }
        if (triple_product(nx2, ny2, nz2) < 0.0) {
            std::cout << "\tWrong face orientation (right-top-front)\n";
            print_cell_info(cell);
            return -1;
        }
    }
    return 0;
}

int check_base_vertices_order(AmrStorage::Item &cell) {
    using zephyr::math::geom::LargeList2D;
    using zephyr::math::geom::LargeList3D;
    using topology::iww;
    using Vector3d;
    using side;

    auto dim = cell[element].dimension;

    auto faces = cell[faces];

    if (dim == 2) {
        LargeList2D vertices;
        for (int i = 0; i < 9; ++i) {
            vertices[i] = (Vector3d &) cell.vertices[i];
        }

        bool bad = false;

        // Индекс пересечения двух граней
        auto cross_face = [](Faces& faces, side side1, side side2) -> int {
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

                if (cross_product(a, b).z < 0.0 ||
                    cross_product(c, d).z < 0.0) {
                    bad = true;
                    break;
                }
            }
            if (bad) {
                break;
            }
        }

        // Пересечения граней по нужным вершинам
        if (cross_face(faces, Side::::LEFT0, Side::BOTTOM0) != iww(0, 0)) {
            bad = true;
        }
        if (cross_face(faces, Side::::LEFT0, Side::TOP) != iww(0, 2) &&
            cross_face(faces, Side::::LEFT1, Side::TOP) != iww(0, 2)) {
            bad = true;
        }
        if (cross_face(faces, Side::::RIGHT, Side::BOTTOM0) != iww(2, 0) &&
            cross_face(faces, Side::::RIGHT, Side::BOTTOM1) != iww(2, 0)) {
            bad = true;
        }
        if (cross_face(faces, Side::::RIGHT0, Side::TOP0) != iww(2, 2) &&
            cross_face(faces, Side::::RIGHT0, Side::TOP1) != iww(2, 2) &&
            cross_face(faces, Side::::RIGHT1, Side::TOP0) != iww(2, 2) &&
            cross_face(faces, Side::::RIGHT1, Side::TOP1) != iww(2, 2)) {
            bad = true;
        }

        if (bad) {
            std::cout << "\tBad arrangement of vertices in cell\n";
            print_cell_info(cell);
            return -1;
        }
    } else {
        LargeList3D vertices;
        for (int i = 0; i < 27; ++i) {
            vertices[i] = (Vector3d &) cell.vertices[i];
        }

        bool bad = false;

        // Индекс пересечения трех граней
        auto cross_face = [](Faces& faces, side side1, side side2, side side3) -> int {
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

                    if (triple_product(a, b, c) < 0.0 ||
                        triple_product(A, B, C) < 0.0) {
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
            print_cell_info(cell);
            return -1;
        }
    }
    return 0;
}

int check_complex_faces(AmrStorage::Item &cell) {
    using Vector3d;
    using faces;
    using vertices;
    using side;

    auto dim = cell[element].dimension;

    if (dim == 2) {
        for (int s = 0; s < 4; ++s) {
            auto f1 = cell.faces[s];
            auto f2 = cell.faces[s + 6];
            if (f2.is_undefined()) continue;

            if (fabs(scalar_product(f1.normal, f2.normal) - 1.0) > 1.0e-2) {
                std::cout << "\tSubfaces are not co-directed (" + side_to_string(side(s)) + " side)\n";
                print_cell_info(cell);
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
                print_cell_info(cell);
                return -1;
            }
        }
    } else {
        for (int s = 0; s < 6; ++s) {
            auto f1 = cell.faces[s];
            auto f2 = cell.faces[s + 6];
            if (f2.is_undefined()) continue;

            auto f3 = cell.faces[s + 12];
            auto f4 = cell.faces[s + 18];
            if (f3.is_undefined() || f4.is_undefined()) {
                std::cout
                        << "\tComplex 3D face (" + side_to_string(side(s)) + " side) has less than 4 subfaces\n";
                print_cell_info(cell);
                return -1;
            }

            double d1 = fabs(scalar_product(f1.normal, f2.normal) - 1.0);
            double d2 = fabs(scalar_product(f1.normal, f3.normal) - 1.0);
            double d3 = fabs(scalar_product(f1.normal, f4.normal) - 1.0);
            if (d1 + d2 + d3 > 1.0e-5) {
                std::cout << "\tSubfaces are not co-directed (" + side_to_string(side(s)) + " side)\n";
                print_cell_info(cell);
                return -1;
            }
        }
    }
    return 0;
}

int check_connectivity(AmrStorage &locals, size_t ic, AmrStorage& aliens, int rank) {
    using Vector3d;

    auto cell = locals[ic];
    if (cell[element].kind == kind::UNDEFINED)
        return 0;

    auto dim = cell[element].dimension;

    // Через обычные грани существуют соседи
    for (int iface = 0; iface < Faces::max_size; ++iface) {
        auto &face = cell.faces[iface];
        if (face.is_undefined()) continue;

        if (face.boundary != Boundary::ORDINARY &&
            face.boundary != Boundary::PERIODIC) {
            // Простая граничная грань
            continue;
        }

        auto neib = locals[0];
        auto neib_rank = face.adjacent.rank;
        auto neib_index = face.adjacent.index;
        auto ghost = face.adjacent.ghost;

        if (neib_rank == rank) {
            // Локальная ячейка
            if (ghost < std::numeric_limits<int>::max()) {
                std::cout << "\tadjacent.ghost != infinity\n";
                print_cell_info(cell);
                return -1;
            }
            if (neib_index >= locals.size()) {
                std::cout << "\tLocal neighbor out of range\n";
                print_cell_info(cell);
                return -1;
            }
            if (neib_index == ic) {
                std::cout << "\tSelf reference\n";
                print_cell_info(cell);
                return -1;
            }
            neib = locals[neib_index];
        }
        else {
            // Удаленная ячейка
            if (ghost >= aliens.size()) {
                std::cout << "\tRemote neighbor out of range\n";
                print_cell_info(cell);
                return -1;
            }

            neib = aliens[ghost];
        }

        if (neib[element].kind == kind::UNDEFINED) {
            std::cout << "\tUndefined neighbor\n";
            print_cell_info(cell);
            return -1;
        }

        Vector3d fc(0.0, 0.0, 0.0);
        for (int i = 0; i < VpF(dim); ++i) {
            fc += (Vector3d &) cell.vertices[face.vertices[i]];
        }
        fc /= VpF(dim);

        // Сосед должен иметь точно такую же грань, но с противоположной нормалью,
        // и ссылаться на текущую ячейку
        int counter = 0;
        for (auto &nface: neib.faces) {
            if (nface.is_undefined()) continue;

            if (nface.boundary != Boundary::ORDINARY &&
                nface.boundary != Boundary::PERIODIC) {
                // простая граничная грань
                continue;
            }

            Vector3d nfc(0.0, 0.0, 0.0);
            for (int i = 0; i < VpF(dim); ++i) {
                nfc += (Vector3d &) neib.vertices[nface.vertices[i]];
            }
            nfc /= VpF(dim);

            if (distance(fc, nfc) > 1.0e-6 * cell[size].value) {
                continue;
            }

            // мы нашли соответствующую грань соседа
            ++counter;

            // нормали противоположны
            if (fabs(scalar_product(face.normal, nface.normal) + 1.0) > 1.0e-6) {
                std::cout << "\tOpposite faces have not opposite normals\n";
                std::cout << "\tCurrent cell:\n";
                print_cell_info(cell);
                std::cout << "\tNeighbor:\n";
                print_cell_info(neib);
                return -1;
            }
            // площади совпадают
            if (fabs(face.area - nface.area) > 1.0e-6 * cell[size].value) {
                std::cout << "\tOpposite faces have different area\n";
                std::cout << "\tCurrent cell:\n";
                print_cell_info(cell);
                std::cout << "\tNeighbor:\n";
                print_cell_info(neib);
                return -1;
            }
            if (nface.adjacent.rank != rank) {
                std::cout << "\tWrong adjacent (rank)\n";
                std::cout << "\tCurrent cell:\n";
                print_cell_info(cell);
                std::cout << "\tNeighbor:\n";
                print_cell_info(neib);
                return -1;
            }
            if (neib_rank == rank) {
                // Локальный сосед
                if (nface.adjacent.ghost < std::numeric_limits<int>::max()) {
                    std::cout << "\tWrong connection (ghost < infinity)\n";
                    std::cout << "\tCurrent cell:\n";
                    print_cell_info(cell);
                    std::cout << "\tNeighbor:\n";
                    print_cell_info(neib);
                    return -1;
                }
            }
            else {
                // Удаленный сосед
                if (nface.adjacent.ghost > std::numeric_limits<int>::max()) {
                    std::cout << "\tWrong connection (huge ghost)\n";
                    std::cout << "\tCurrent cell:\n";
                    print_cell_info(cell);
                    std::cout << "\tNeighbor:\n";
                    print_cell_info(neib);
                    return -1;
                }

            }
        }
        if (counter < 1) {
            std::cout << "\tHas no neighbor across ordinary face (" << side_to_string(iface/4) << ")\n";
            print_cell_info(cell);
            std::cout << "\tNeighbor:\n";
            print_cell_info(locals[cell.faces[iface].adjacent.index]);
            return -1;
        }
        if (counter > 1) {
            std::cout << "\tMore than one neighbor across ordinary face\n";
            print_cell_info(cell);
            return -1;
        }
    }

    return 0;
}

int check_base_mesh(AmrStorage &locals, AmrStorage &aliens, int rank) {
    using Vector3d;
    using coords;
    using vertices;
    using faces;
    using element;
    using side;

    if (locals.empty()) {
        return 0;
    }

    test_storage(locals);

    size_t n_cells = locals.size();
    if (n_cells < 1) {
        std::cout << "\tEmpty storage\n";
        return -1;
    }

    auto dim = locals[0][element].dimension;
    auto kind = locals[0][element].kind;

    if (dim != 2 && dim != 3) {
        std::cout << "\tDimension is not 2 or 3\n";
        return -1;
    }

    if (kind != kind::EULER) {
        std::cout << "\tKind is not EULER\n";
        return -1;
    }

    int res = 0;
    for (size_t ic = 0; ic < locals.size(); ++ic) {
        auto cell = locals[ic];

        // Размерность постоянна
        if (cell[element].dimension != dim) {
            std::cout << "\tVarious dimensions of elements\n";
            return -1;
        }

        // Ячейки одно типа
        if (cell[element].kind != kind) {
            std::cout << "\tVarious kinds of elements\n";
            return -1;
        }

        // Число граней
        for (int i = 0; i < FpC(dim); ++i) {
            if (cell.faces[i].is_undefined()) {
                std::cout << "\tCell has no one of main faces\n";
                print_cell_info(cell);
                return -1;
            }
        }
        if (cell[faces].size() > FpC(dim)) {
            std::cout << "\tCell has too much faces\n";
            print_cell_info(cell);
            return -1;
        }

        // Число вершин ???

        // Правильное задание геометрии
        res = check_geometry(cell);
        if (res < 0) return res;

        // Грани правльно ориентированы
        res = check_base_face_orientation(cell);
        if (res < 0) return res;

        // Порядок основных вершин
        res = check_base_vertices_order(cell);
        if (res < 0) return res;

        // Проверка смежности
        res = check_connectivity(locals, ic, aliens, rank);
        if (res < 0) return res;
    }

    return 0;
}

int check_refined_mesh(zephyr::data::AmrStorage &locals, AmrStorage &aliens, int rank) {
    using Vector3d;
    using coords;
    using vertices;
    using faces;
    using element;
    using side;

    if (locals.empty()) {
        return 0;
    }

    test_storage(locals);

    size_t n_cells = locals.size();
    if (n_cells < 1) {
        std::cout << "\tEmpty storage\n";
        return -1;
    }

    auto dim = locals[0][element].dimension;
    auto kind = locals[0][element].kind;

    if (dim != 2 && dim != 3) {
        std::cout << "\tDimension is not 2 or 3\n";
        return -1;
    }

    if (kind != kind::EULER) {
        std::cout << "\tKind is not EULER\n";
        return -1;
    }

    int res = 0;
    for (size_t ic = 0; ic < locals.size(); ++ic) {
        auto cell = locals[ic];

        if (cell[element].is_undefined()) {
            continue;
        }

        // Размерность постоянна
        if (cell[element].dimension != dim) {
            std::cout << "\tVarious dimensions of elements\n";
            std::cout << "\t\tdimension: " << dim << "\n";
            std::cout << "\t\tcell.dimension: " << cell[element].dimension << "\n";
            return -1;
        }

        // Ячейки одно типа
        if (cell[element].kind != kind && cell[element].kind != kind::UNDEFINED) {
            std::cout << "\tVarious kinds of elements\n";
            return -1;
        }

        // Число граней
        for (int i = 0; i < FpC(dim); ++i) {
            if (cell.faces[i].is_undefined()) {
                std::cout << "\tCell has no one of main faces\n";
                print_cell_info(cell);
                return -1;
            }
        }

        // Вершины дублируются
        for (int i = 0; i < cell[vertices].size(); ++i) {
            for (int j = i + 1; j < cell[vertices].size(); ++j) {
                double dist = distance(
                        (Vector3d &) cell.vertices[i],
                        (Vector3d &) cell.vertices[j]);
                if (dist < 1.0e-5 * cell[size]) {
                    std::cout << "\tIdentical vertices\n";
                    print_cell_info(cell);
                    return -1;
                }
            }
        }

        // Правильное задание геометрии
        res = check_geometry(cell);
        if (res < 0) return res;

        // Грани правльно ориентированы
        res = check_base_face_orientation(cell);
        if (res < 0) return res;

        // Порядок основных вершин
        res = check_base_vertices_order(cell);
        if (res < 0) return res;

        // Проверка сложных граней
        res = check_complex_faces(cell);
        if (res < 0) return res;

        // Проверка смежности
        res = check_connectivity(locals, ic, aliens, rank);
        if (res < 0) return res;
    }

    return 0;
}

} // namespace amr2
} // namespace mesh
} // namespace zephyr
 */