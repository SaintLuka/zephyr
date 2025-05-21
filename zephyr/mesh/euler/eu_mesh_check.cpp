#include <zephyr/utils/mpi.h>

#include <zephyr/mesh/primitives/side.h>
#include <zephyr/mesh/primitives/amr_cell.h>
#include <zephyr/mesh/primitives/bfaces.h>
#include <zephyr/mesh/primitives/decomposition.h>
#include <zephyr/mesh/euler/eu_mesh.h>

namespace zephyr::mesh {

using utils::mpi;
using namespace geom;

int check_connectivity(AmrStorage &locals, int ic, AmrStorage& aliens) {
    auto& cell = locals[ic];
    if (cell.is_undefined()) {
        return 0;
    }

    auto dim = cell.dim;

    // Через обычные грани существуют соседи
    for (int iface = 0; iface < BFaces::max_count; ++iface) {
        auto &face = cell.faces[iface];
        if (face.is_undefined()) {
            // Проверим обнуление параметров
            if (face.adjacent.rank >= 0) {
                std::cout << "\tUndefined face, adjacent.rank >= 0\n";
                cell.print_info();
                return -1;
            }
            if (face.adjacent.alien >= 0) {
                std::cout << "\tUndefined face, adjacent.alien >= 0\n";
                cell.print_info();
                return -1;
            }
            if (face.adjacent.index >= 0) {
                std::cout << "\tUndefined face, adjacent.index >= 0\n";
                cell.print_info();
                return -1;
            }
            continue;
        }

        // Простая граничная грань
        if (face.boundary != Boundary::ORDINARY &&
            face.boundary != Boundary::PERIODIC) {

            // Грань должна ссылаться на искомую ячейку
            if (face.adjacent.rank != cell.rank ||
                face.adjacent.index != cell.index ||
                face.adjacent.alien >= 0) {
                std::cout << "\tBoundary face should point to origin cell\n";
                cell.print_info();
                return -1;
            }

            continue;
        }

        AmrCell* neib_ptr = nullptr;
        auto& adj = face.adjacent;

        if (adj.rank < 0) {
            std::cout << "\tadjacent.rank < 0\n";
            cell.print_info();
            return -1;
        }

        if (adj.rank == mpi::rank()) {
            // Локальная ячейка
            if (adj.alien >= 0) {
                std::cout << "\tadjacent.alien >= 0\n";
                cell.print_info();
                return -1;
            }
            if (adj.index < 0 || adj.index >= locals.size()) {
                std::cout << "\tLocal neighbor out of range\n";
                cell.print_info();
                return -1;
            }
            if (adj.index == ic) {
                std::cout << "\tSelf reference\n";
                cell.print_info();
                return -1;
            }
            neib_ptr = &locals[adj.index];
        }
        else {
            // Удаленная ячейка
            if (adj.index < 0) {
                std::cout << "\tRemote neighbor wrong index\n";
                cell.print_info();
                return -1;
            }
            if (adj.alien < 0 || adj.alien >= aliens.size()) {
                std::cout << "\tRemote neighbor out of range\n";
                cell.print_info();
                return -1;
            }

            neib_ptr = &aliens[adj.alien];
        }

        if (!neib_ptr || neib_ptr->is_undefined()) {
            std::cout << "\tUndefined neighbor\n";
            cell.print_info();
            return -1;
        }

        AmrCell& neib = *neib_ptr;

        Vector3d fc(0.0, 0.0, 0.0);
        for (int i = 0; i < VpF(dim); ++i) {
            fc += cell.vertices[face.vertices[i]];
        }
        fc /= VpF(dim);

        // Сосед должен иметь точно такую же грань, но с противоположной нормалью,
        // и ссылаться на текущую ячейку
        int counter = 0;
        for (auto &nface: neib.faces) {
            if (nface.is_undefined()) {
                continue;
            }

            if (nface.boundary != Boundary::ORDINARY &&
                nface.boundary != Boundary::PERIODIC) {
                // простая граничная грань
                continue;
            }

            Vector3d nfc(0.0, 0.0, 0.0);
            for (int i = 0; i < VpF(dim); ++i) {
                nfc += neib.vertices[nface.vertices[i]];
            }
            nfc /= VpF(dim);

            if ((fc - nfc).norm() > 1.0e-6 * cell.linear_size()) {
                continue;
            }

            // мы нашли соответствующую грань соседа
            ++counter;

            // нормали противоположны
            if (std::abs(face.normal.dot(nface.normal) + 1.0) > 1.0e-6) {
                std::cout << "\tOpposite faces have not opposite normals\n";
                std::cout << "\tCurrent cell:\n";
                cell.print_info();
                std::cout << "\tNeighbor:\n";
                neib.print_info();
                return -1;
            }
            // площади совпадают
            if (std::abs(face.area - nface.area) > 1.0e-6 * cell.linear_size()) {
                std::cout << "\tOpposite faces have different area\n";
                std::cout << "\tCurrent cell:\n";
                cell.print_info();
                std::cout << "\tNeighbor:\n";
                neib.print_info();
                return -1;
            }
            // Указывает на исходную ячейку
            if (nface.adjacent.index != ic) {
                std::cout << "\tWrong connection (index != ic)\n";
                std::cout << "\tCurrent cell:\n";
                cell.print_info();
                std::cout << "\tNeighbor:\n";
                neib.print_info();
                return -1;
            }
            // Ранг смежной (исходная) больше нуля
            if (nface.adjacent.rank < 0) {
                std::cout << "\tWrong adjacent (rank < 0)\n";
                std::cout << "\tCurrent cell:\n";
                cell.print_info();
                std::cout << "\tNeighbor:\n";
                neib.print_info();
                return -1;
            }
            // Ранг смежной (исходная) равен рангу процесса
            if (nface.adjacent.rank != mpi::rank()) {
                std::cout << "\tWrong adjacent (rank)\n";
                std::cout << "\tCurrent cell:\n";
                cell.print_info();
                std::cout << "\tNeighbor:\n";
                neib.print_info();
                return -1;
            }
            if (adj.rank == mpi::rank()) {
                // Локальный сосед
                if (nface.adjacent.alien >= 0) {
                    std::cout << "\tWrong connection (alien >= 0)\n";
                    std::cout << "\tCurrent cell:\n";
                    cell.print_info();
                    std::cout << "\tNeighbor:\n";
                    neib.print_info();
                    return -1;
                }
            }
            else {
                // Удаленный сосед
                if (nface.adjacent.alien < 0) {
                    // Проверка не актуальна сразу после построения списка соседей,
                    // face.alien для alien ячеек не получены, да и не важны.
                    /*
                    std::cout << "\tWrong connection\n";
                    std::cout << "\tCurrent cell:\n";
                    cell.print_info();
                    std::cout << "\tNeighbor:\n";
                    neib.print_info();
                    return -1;
                     */
                }
            }
        }
        if (counter < 1) {
            std::cout << "\tHas no neighbor across ordinary " << side_to_string(iface, dim) << ")\n";
            cell.print_info();
            std::cout << "\tNeighbor:\n";
            locals[cell.faces[iface].adjacent.index].print_info();
            return -1;
        }
        if (counter > 1) {
            std::cout << "\tMore than one neighbor across ordinary face\n";
            cell.print_info();
            return -1;
        }
    }

    return 0;
}

int EuMesh::check_base() {
    if (m_locals.empty()) {
        if (mpi::single()) {
            std::cout << "\tEmpty storage\n";
            return -1;
        } else {
            return 0;
        }
    }

    auto dim = m_locals[0].dim;

    if (dim != 2 && dim != 3) {
        std::cout << "\tDimension is not 2 or 3\n";
        return -1;
    }

    int res = 0;
    for (int ic = 0; ic < m_locals.size(); ++ic) {
        AmrCell& cell = m_locals[ic];

        if (cell.index < 0 || cell.index != ic) {
            std::cout << "\tWrong cell index\n";
            return -1;
        }

        if (cell.rank < 0 || cell.rank != mpi::rank()) {
            std::cout << "\tWrong cell rank\n";
            return -1;
        }

        // Размерность постоянна
        if (cell.dim != dim) {
            std::cout << "\tVarious dimensions of elements\n";
            return -1;
        }

        // Число граней
        for (int i = 0; i < FpC(dim); ++i) {
            if (cell.faces[i].is_undefined()) {
                std::cout << "\tCell has no one of main faces\n";
                cell.print_info();
                return -1;
            }
        }
        if (cell.faces.count() > FpC(dim)) {
            std::cout << "\tCell has too much faces\n";
            cell.print_info();
            return -1;
        }

        // Число вершин ???

        // Правильное задание геометрии
        res = cell.check_geometry();
        if (res < 0) return res;

        // Грани правльно ориентированы
        res = cell.check_base_face_orientation();
        if (res < 0) return res;

        // Порядок основных вершин
        res = cell.check_base_vertices_order();
        if (res < 0) return res;

        // Проверка смежности
        res = check_connectivity(m_locals, ic, m_aliens);
        if (res < 0) return res;
    }

    return 0;
}

int EuMesh::check_refined() {
    if (m_locals.empty()) {
        if (mpi::single()) {
            std::cout << "\tEmpty storage\n";
            return -1;
        } else {
            return 0;
        }
    }

    auto dim = m_locals[0].dim;

    if (dim != 2 && dim != 3) {
        std::cout << "\tDimension is not 2 or 3\n";
        return -1;
    }

    int res = 0;
    for (int ic = 0; ic < m_locals.size(); ++ic) {
        AmrCell& cell = m_locals[ic];

        if (cell.is_undefined()) {
            std::cout << "\tUndefined cell\n";
            return -1;
        }

        if (cell.index < 0 || cell.index != ic) {
            std::cout << "\tWrong cell index\n";
            return -1;
        }

        if (cell.rank < 0 || cell.rank != mpi::rank()) {
            std::cout << "\tWrong cell rank\n";
            return -1;
        }

        // Размерность постоянна
        if (cell.dim != dim) {
            std::cout << "\tVarious dimensions of elements\n";
            std::cout << "\t\tdimension: " << dim << "\n";
            std::cout << "\t\tcell.dimension: " << cell.dim << "\n";
            return -1;
        }

        // Число граней
        for (int i = 0; i < FpC(dim); ++i) {
            if (cell.faces[i].is_undefined()) {
                std::cout << "\tCell has no one of main faces\n";
                cell.print_info();
                return -1;
            }
        }

        // Вершины дублируются
        for (int i = 0; i < std::pow(3, dim); ++i) {
            for (int j = i + 1; j < std::pow(3, dim); ++j) {
                double dist = (cell.vertices[i] - cell.vertices[j]).norm();
                if (dist < 1.0e-5 * cell.linear_size()) {
                    std::cout << "\tIdentical vertices\n";
                    cell.print_info();
                    return -1;
                }
            }
        }

        // Правильное задание геометрии
        res = cell.check_geometry();
        if (res < 0) return res;

        // Грани правльно ориентированы
        res = cell.check_base_face_orientation();
        if (res < 0) return res;

        // Порядок основных вершин
        res = cell.check_base_vertices_order();
        if (res < 0) return res;

        // Проверка сложных граней
        res = cell.check_complex_faces();
        if (res < 0) return res;

        // Проверка смежности
        res = check_connectivity(m_locals, ic, m_aliens);
        if (res < 0) return res;
    }

    return 0;
}

} // namespace zephyr::mesh