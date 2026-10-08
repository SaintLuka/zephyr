#include <zephyr/mesh/amr/rotations.h>
#include <zephyr/geom/indexing.h>

using namespace zephyr::mesh::amr;
using namespace zephyr::geom::indexing;
using zephyr::mesh::SqMap;

template <int dim>
SqMap<dim> unit_cell() {
    if constexpr (dim == 2) {
        std::array verts_2D = {
            Vector3d{-0.5, -0.5, 0.0},
            Vector3d{+0.5, -0.5, 0.0},
            Vector3d{-0.5, +0.5, 0.0},
            Vector3d{+0.5, +0.5, 0.0}
        };

        return SqQuad(verts_2D[0], verts_2D[1], verts_2D[2], verts_2D[3]);
    }
    else {
        std::array verts_3D = {
            Vector3d{-0.5, -0.5, -0.5},
            Vector3d{+0.5, -0.5, -0.5},
            Vector3d{-0.5, +0.5, -0.5},
            Vector3d{+0.5, +0.5, -0.5},
            Vector3d{-0.5, -0.5, +0.5},
            Vector3d{+0.5, -0.5, +0.5},
            Vector3d{-0.5, +0.5, +0.5},
            Vector3d{+0.5, +0.5, +0.5}
        };

        return SqCube(verts_3D[0], verts_3D[1], verts_3D[2], verts_3D[3],
                      verts_3D[4], verts_3D[5], verts_3D[6], verts_3D[7]);
    }

}

template <int dim>
void role_in_child_gen() {
    int n_children_at_node[VpC(dim)];
    int children_at_node[VpC(dim)][CpC(dim)];
    int role_in_child[VpC(dim)][CpC(dim)];

    auto cell = unit_cell<dim>();
    auto children = cell.children();

    for (int role = 0; role < VpC(dim); ++role) {
        Vector3d v = cell[role];

        int count = 0;
        for (int ich = 0; ich < CpC(dim); ++ich) {
            const auto& child = children[ich];
            int res = -1;
            for (int r = 0; r < VpC(dim); ++r) {
                if ((v - child[r]).norm() < 1.0e-5) {
                    res = r;
                }
            }
            role_in_child[role][ich] = res;
            if (res >= 0) {
                children_at_node[role][count] = ich;
                ++count;
            }
        }
        for (int i = count; i < CpC(dim); ++i) {
            children_at_node[role][i] = -1;
        }
        n_children_at_node[role] = count;
    }

    std::cout << "/// @brief Число дочерних ячеек при узле [role]\n";
    if constexpr (dim == 2) {
        std::cout << "static constexpr int n_children_at_node_2D[9] = {\n";
    }
    else {
        std::cout << "static constexpr int n_children_at_node_3D[27] = {\n";
    }
    for (int role = 0; role < VpC(dim); ++role) {
        if (role % 3 == 0) std::cout << "    ";
        std::cout << n_children_at_node[role] << ", ";
        if (role % 3 == 2) std::cout << "\n";
    }
    std::cout << "};\n\n";

    std::cout << "/// @brief Индексы дочерних ячеек около узла [role][child]\n";
    if constexpr (dim == 2) {
        std::cout << "static constexpr int children_at_node_2D[9][4] = {\n";
    }
    else {
        std::cout << "static constexpr int children_at_node_3D[27][8] = {\n";
    }
    for (int role = 0; role < VpC(dim); ++role) {
        std::cout << "    {";
        for (int ich = 0; ich < CpC(dim) - 1; ++ich) {
            std::cout << std::setw(2) << children_at_node[role][ich] << ", ";
        }
        std::cout << std::setw(2) << children_at_node[role][CpC(dim) - 1] << "},\n";
    }
    std::cout << "};\n\n";

    std::cout << "/// @brief Роль узла родительской ячейки в дочерней [role][child]\n";
    if constexpr (dim == 2) {
        std::cout << "static constexpr int role_in_child_2D[9][4] = {\n";
    }
    else {
        std::cout << "static constexpr int role_in_child_3D[27][8] = {\n";
    }
    for (int role = 0; role < VpC(dim); ++role) {
        std::cout << "    {";
        for (int ich = 0; ich < CpC(dim) - 1; ++ich) {
            std::cout << std::setw(2) << role_in_child[role][ich] << ", ";
        }
        std::cout << std::setw(2) << role_in_child[role][CpC(dim) - 1] << "},\n";
    }
    std::cout << "};\n\n";
}

int main() {
    role_in_child_gen<2>();
    role_in_child_gen<3>();

    //children_at_node_gen<2>();
    //children_at_node_gen<3>();

    return 0;
}