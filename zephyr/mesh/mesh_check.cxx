#include <zephyr/utils/mpi.h>
#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/mesh.h>

namespace zephyr { namespace mesh {

using utils::mpi;

int check_geometry(const Cell& cell) {
    /*
    auto dim = cell.dim;

    for (auto &face: cell.faces) {
        if (face.is_undefined()) continue;

        Vector3d fc(0.0, 0.0, 0.0);
        for (int iv = 0; iv < amr::VpF(dim); ++iv) {
            fc += cell.vertices[face.vertices[iv]];
        }
        fc /= amr::VpF(dim);

        // Нормаль внешняя
        if (face.normal.dot(fc - cell.coords) < 0.0) {
            std::cout << "\tWrong normal direction (inside cell)\n";
            cell.print_info();
            return -1;
        }

        // Вершины грани перечислены в правильном порядке
        if (dim > 2) {
            Vector3d v0 = cell.vertices[face.vertices[0]];
            Vector3d v1 = cell.vertices[face.vertices[1]];
            Vector3d v2 = cell.vertices[face.vertices[2]];
            Vector3d v3 = cell.vertices[face.vertices[3]];

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
    */
}

int check_base_face_orientation(Storage::Item &cell) {
    /*

    if (cell[element].dimension == 2) {
        Vector3d nx1 = (Vector3d &) cell.geom().faces[Side::LEFT].normal;
        Vector3d nx2 = (Vector3d &) cell.geom().faces[Side::RIGHT].normal;
        Vector3d ny1 = (Vector3d &) cell.geom().faces[Side::BOTTOM].normal;
        Vector3d ny2 = (Vector3d &) cell.geom().faces[Side::TOP].normal;

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
        Vector3d nx1 = (Vector3d &) cell.geom().faces[Side::LEFT].normal;
        Vector3d nx2 = (Vector3d &) cell.geom().faces[Side::RIGHT].normal;
        Vector3d ny1 = (Vector3d &) cell.geom().faces[Side::BOTTOM].normal;
        Vector3d ny2 = (Vector3d &) cell.geom().faces[Side::TOP].normal;
        Vector3d nz1 = (Vector3d &) cell.geom().faces[Side::BACK].normal;
        Vector3d nz2 = (Vector3d &) cell.geom().faces[Side::FRONT].normal;

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
     */
    return 0;
}

int check_base_vertices_order(Storage::Item &cell) {
    /*
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
            vertices[i] = (Vector3d &) cell.geom().vertices[i];
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
            vertices[i] = (Vector3d &) cell.geom().vertices[i];
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
     */

    return 0;
}

int check_complex_faces(Storage::Item &cell) {
    /*
    using Vector3d;
    using faces;
    using vertices;
    using side;

    auto dim = cell[element].dimension;

    if (dim == 2) {
        for (int s = 0; s < 4; ++s) {
            auto f1 = cell.geom().faces[s];
            auto f2 = cell.geom().faces[s + 6];
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
            auto f1 = cell.geom().faces[s];
            auto f2 = cell.geom().faces[s + 6];
            if (f2.is_undefined()) continue;

            auto f3 = cell.geom().faces[s + 12];
            auto f4 = cell.geom().faces[s + 18];
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

     */
    return 0;
}

int check_connectivity(Storage &locals, int ic, Storage& aliens) {
    /*
    using Vector3d;

    auto cell = locals[ic];
    if (cell[element].kind == kind::UNDEFINED)
        return 0;

    auto dim = cell[element].dimension;

    // Через обычные грани существуют соседи
    for (int iface = 0; iface < Faces::max_size; ++iface) {
        auto &face = cell.geom().faces[iface];
        if (face.is_undefined()) continue;

        if (face.boundary != FaceFlag::ORDINARY &&
            face.boundary != FaceFlag::PERIODIC) {
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
            fc += (Vector3d &) cell.geom().vertices[face.vertices[i]];
        }
        fc /= VpF(dim);

        // Сосед должен иметь точно такую же грань, но с противоположной нормалью,
        // и ссылаться на текущую ячейку
        int counter = 0;
        for (auto &nface: neib.geom().faces) {
            if (nface.is_undefined()) continue;

            if (nface.boundary != FaceFlag::ORDINARY &&
                nface.boundary != FaceFlag::PERIODIC) {
                // простая граничная грань
                continue;
            }

            Vector3d nfc(0.0, 0.0, 0.0);
            for (int i = 0; i < VpF(dim); ++i) {
                nfc += (Vector3d &) neib.geom().vertices[nface.vertices[i]];
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
            print_cell_info(locals[cell.geom().faces[iface].adjacent.index]);
            return -1;
        }
        if (counter > 1) {
            std::cout << "\tMore than one neighbor across ordinary face\n";
            print_cell_info(cell);
            return -1;
        }
    }
     */

    return 0;
}

int Mesh::check_base() {
    if (m_locals.empty()) {
        if (mpi::is_single()) {
            std::cout << "\tEmpty storage\n";
            return -1;
        } else {
            return 0;
        }
    }

    auto dim = m_locals[0].dim();

    if (dim != 2 && dim != 3) {
        std::cout << "\tDimension is not 2 or 3\n";
        return -1;
    }

    int res = 0;
    for (int ic = 0; ic < m_locals.size(); ++ic) {
        auto cell = m_locals[ic];

        // Размерность постоянна
        if (cell.dim() != dim) {
            std::cout << "\tVarious dimensions of elements\n";
            return -1;
        }

        // Число граней
        for (int i = 0; i < amr::FpC(dim); ++i) {
            if (cell.faces(i).is_undefined()) {
                std::cout << "\tCell has no one of main faces\n";
                cell.print_info();
                return -1;
            }
        }
        if (cell.faces().size() > amr::FpC(dim)) {
            std::cout << "\tCell has too much faces\n";
            cell.print_info();
            return -1;
        }

        // Число вершин ???

        // Правильное задание геометрии
        res = check_geometry(cell.geom());
        if (res < 0) return res;

        // Грани правльно ориентированы
        res = check_base_face_orientation(cell);
        if (res < 0) return res;

        // Порядок основных вершин
        res = check_base_vertices_order(cell);
        if (res < 0) return res;

        // Проверка смежности
        res = check_connectivity(m_locals, ic, m_aliens);
        if (res < 0) return res;
    }

    return 0;
}

int Mesh::check_refined() {
    if (m_locals.empty()) {
        if (mpi::is_single()) {
            std::cout << "\tEmpty storage\n";
            return -1;
        } else {
            return 0;
        }
    }

    auto dim = m_locals[0].dim();

    if (dim != 2 && dim != 3) {
        std::cout << "\tDimension is not 2 or 3\n";
        return -1;
    }

    int res = 0;
    for (int ic = 0; ic < m_locals.size(); ++ic) {
        auto cell = m_locals[ic];

        if (cell.is_undefined()) {
            continue;
        }

        // Размерность постоянна
        if (cell.dim() != dim) {
            std::cout << "\tVarious dimensions of elements\n";
            std::cout << "\t\tdimension: " << dim << "\n";
            std::cout << "\t\tcell.dimension: " << cell.dim() << "\n";
            return -1;
        }

        // Число граней
        for (int i = 0; i < amr::FpC(dim); ++i) {
            if (cell.faces(i).is_undefined()) {
                std::cout << "\tCell has no one of main faces\n";
                cell.print_info();
                return -1;
            }
        }

        // Вершины дублируются
        for (int i = 0; i < cell.geom().vertices.size(); ++i) {
            for (int j = i + 1; j < cell.geom().vertices.size(); ++j) {
                double dist = (cell.vertices(i) - cell.vertices(j)).norm();
                if (dist < 1.0e-5 * cell.size()) {
                    std::cout << "\tIdentical vertices\n";
                    cell.print_info();
                    return -1;
                }
            }
        }

        // Правильное задание геометрии
        res = check_geometry(cell.geom());
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
        res = check_connectivity(m_locals, ic, m_aliens);
        if (res < 0) return res;
    }

    return 0;
}

} // mesh
} // zephyr