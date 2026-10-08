// Не устанавливается при установке zephyr, детали алгоритмов и комментарии
// к функциям предназначены для разработчиков.
#pragma once

#include <zephyr/mesh/amr/common.h>
#include <zephyr/mesh/amr/statistics.h>

namespace zephyr::mesh::amr {

/// @brief Число дочерних ячеек при узле [role]
static constexpr int n_children_at_node_2D[9] = {
    1, 2, 1,
    2, 4, 2,
    1, 2, 1,
};

/// @brief Индексы дочерних ячеек около узла [role][child]
static constexpr int children_at_node_2D[9][4] = {
    { 0, -1, -1, -1},
    { 0,  1, -1, -1},
    { 1, -1, -1, -1},
    { 0,  2, -1, -1},
    { 0,  1,  2,  3},
    { 1,  3, -1, -1},
    { 2, -1, -1, -1},
    { 2,  3, -1, -1},
    { 3, -1, -1, -1},
};

/// @brief Роль узла родительской ячейки в дочерней [role][child]
static constexpr int role_in_child_2D[9][4] = {
    { 0, -1, -1, -1},
    { 2,  0, -1, -1},
    {-1,  2, -1, -1},
    { 6, -1,  0, -1},
    { 8,  6,  2,  0},
    {-1,  8, -1,  2},
    {-1, -1,  6, -1},
    {-1, -1,  8,  6},
    {-1, -1, -1,  8},
};

/// @brief Число дочерних ячеек при узле [role]
static constexpr int n_children_at_node_3D[27] = {
    1, 2, 1,
    2, 4, 2,
    1, 2, 1,
    2, 4, 2,
    4, 8, 4,
    2, 4, 2,
    1, 2, 1,
    2, 4, 2,
    1, 2, 1,
};

/// @brief Индексы дочерних ячеек около узла [role][child]
static constexpr int children_at_node_3D[27][8] = {
    { 0, -1, -1, -1, -1, -1, -1, -1},
    { 0,  1, -1, -1, -1, -1, -1, -1},
    { 1, -1, -1, -1, -1, -1, -1, -1},
    { 0,  2, -1, -1, -1, -1, -1, -1},
    { 0,  1,  2,  3, -1, -1, -1, -1},
    { 1,  3, -1, -1, -1, -1, -1, -1},
    { 2, -1, -1, -1, -1, -1, -1, -1},
    { 2,  3, -1, -1, -1, -1, -1, -1},
    { 3, -1, -1, -1, -1, -1, -1, -1},
    { 0,  4, -1, -1, -1, -1, -1, -1},
    { 0,  1,  4,  5, -1, -1, -1, -1},
    { 1,  5, -1, -1, -1, -1, -1, -1},
    { 0,  2,  4,  6, -1, -1, -1, -1},
    { 0,  1,  2,  3,  4,  5,  6,  7},
    { 1,  3,  5,  7, -1, -1, -1, -1},
    { 2,  6, -1, -1, -1, -1, -1, -1},
    { 2,  3,  6,  7, -1, -1, -1, -1},
    { 3,  7, -1, -1, -1, -1, -1, -1},
    { 4, -1, -1, -1, -1, -1, -1, -1},
    { 4,  5, -1, -1, -1, -1, -1, -1},
    { 5, -1, -1, -1, -1, -1, -1, -1},
    { 4,  6, -1, -1, -1, -1, -1, -1},
    { 4,  5,  6,  7, -1, -1, -1, -1},
    { 5,  7, -1, -1, -1, -1, -1, -1},
    { 6, -1, -1, -1, -1, -1, -1, -1},
    { 6,  7, -1, -1, -1, -1, -1, -1},
    { 7, -1, -1, -1, -1, -1, -1, -1},
};

/// @brief Роль узла родительской ячейки в дочерней [role][child]
static constexpr int role_in_child_3D[27][8] = {
    { 0, -1, -1, -1, -1, -1, -1, -1},
    { 2,  0, -1, -1, -1, -1, -1, -1},
    {-1,  2, -1, -1, -1, -1, -1, -1},
    { 6, -1,  0, -1, -1, -1, -1, -1},
    { 8,  6,  2,  0, -1, -1, -1, -1},
    {-1,  8, -1,  2, -1, -1, -1, -1},
    {-1, -1,  6, -1, -1, -1, -1, -1},
    {-1, -1,  8,  6, -1, -1, -1, -1},
    {-1, -1, -1,  8, -1, -1, -1, -1},
    {18, -1, -1, -1,  0, -1, -1, -1},
    {20, 18, -1, -1,  2,  0, -1, -1},
    {-1, 20, -1, -1, -1,  2, -1, -1},
    {24, -1, 18, -1,  6, -1,  0, -1},
    {26, 24, 20, 18,  8,  6,  2,  0},
    {-1, 26, -1, 20, -1,  8, -1,  2},
    {-1, -1, 24, -1, -1, -1,  6, -1},
    {-1, -1, 26, 24, -1, -1,  8,  6},
    {-1, -1, -1, 26, -1, -1, -1,  8},
    {-1, -1, -1, -1, 18, -1, -1, -1},
    {-1, -1, -1, -1, 20, 18, -1, -1},
    {-1, -1, -1, -1, -1, 20, -1, -1},
    {-1, -1, -1, -1, 24, -1, 18, -1},
    {-1, -1, -1, -1, 26, 24, 20, 18},
    {-1, -1, -1, -1, -1, 26, -1, 20},
    {-1, -1, -1, -1, -1, -1, 24, -1},
    {-1, -1, -1, -1, -1, -1, 26, 24},
    {-1, -1, -1, -1, -1, -1, -1, 26},
};

template <int dim>
std::span<const int> children_at_node(int role) {
    if constexpr (dim == 2) {
        return std::span(children_at_node_2D[role], n_children_at_node_2D[role]);
    }
    else {
        return std::span(children_at_node_3D[role], n_children_at_node_3D[role]);
    }
}

template<int dim>
constexpr role_t role_in_child(role_t role, int child) {
    if constexpr (dim == 2) {
        return role_in_child_2D[role][child];
    }
    return role_in_child_3D[role][child];
}

template<int dim>
bool is_center_node(role_t role) {
    static constexpr std::array res_2D = {
        false, false, false,
        false, true , false,
        false, false, false,
    };
    static constexpr std::array res_3D = {
        false, false, false,
        false, false, false,
        false, false, false,

        false, false, false,
        false, true , false,
        false, false, false,

        false, false, false,
        false, false, false,
        false, false, false,
    };
    if constexpr (dim == 2)
        return res_2D[role];
    return res_3D[role];
}

template<int dim>
bool is_face_node(role_t role) {
    static constexpr std::array res_2D = {
        false, true , false,
        true , false, true ,
        false, true , false,
    };
    static constexpr std::array res_3D = {
        false, false, false,
        false, true , false,
        false, false, false,

        false, true , false,
        true , false, true ,
        false, true , false,

        false, false, false,
        false, true , false,
        false, false, false,
    };
    if constexpr (dim == 2)
        return res_2D[role];
    return res_3D[role];
}

template<int dim>
bool is_edge_node(role_t role) {
    static constexpr std::array res_3D = {
        false, true , false,
        true , false, true ,
        false, true , false,

        true , false, true ,
        false, false, false,
        true , false, true ,

        false, true , false,
        true , false, true ,
        false, true , false,
    };
    if constexpr (dim == 2)
        return false;
    return res_3D[role];
}

template<int dim>
bool is_corner_node(role_t role) {
    static constexpr std::array res_2D = {
        true , false, true ,
        false, false, false,
        true , false, true ,
    };
    static constexpr std::array res_3D = {
        true , false, true ,
        false, false, false,
        true , false, true ,

        false, false, false,
        false, false, false,
        false, false, false,

        true , false, true ,
        false, false, false,
        true , false, true ,
    };
    if constexpr (dim == 2)
        return res_2D[role];
    return res_3D[role];
}

template<int dim>
void update_incident(const RawCells& cells, RawNodes& nodes,
    const Statistics& count, const std::vector<index_t>& split_indices) {

    constexpr int inc_count = RawIncident::max_incident_amr(dim);

    index_t first_new_node = nodes.n_nodes();
    index_t nodes_per_split = dim == 2 ? (5*5 - 3*3) : (5*5*5 - 3*3*3);
    nodes.resize_amr(nodes.n_nodes() + split_indices.size() * nodes_per_split, dim);

    for (index_t inode = 0; inode < nodes.n_nodes(); ++inode) {
        index_t offset = nodes.incident.offsets[inode];

        int new_count = 0;
        std::array<role_t,  inc_count> new_roles; new_roles.fill(-1);
        std::array<int,     inc_count> new_ranks; new_ranks.fill(mpi::rank());
        std::array<index_t, inc_count> new_index; new_index.fill(-1);
        std::array<index_t, inc_count> new_ghost; new_ghost.fill(-1);

        for (int j = 0; j < inc_count; ++j) {
            int role = nodes.incident.role[offset + j];
            if (role < 0) break;

            index_t index = nodes.incident.index[offset + j];

            // Инцидентная ячейка ничего не делает, но может переместиться
            if (cells.flag[index] == 0) {
                new_roles[new_count] = role;
                new_ranks[new_count] = 0;
                new_index[new_count] = cells.next[index];
                new_ghost[new_count] = -1;
                ++new_count;
                continue;
            }

            // Инцидентная ячейка бьется, надо выяснить роль узла в каждой дочерней
            // ячейке, и все дочерние ячейки при данном узле записать в инцидентные
            if (cells.flag[index] > 0) {
                for (int ich: children_at_node<dim>(role)) {
                    new_roles[new_count] = role_in_child<dim>(role, ich);
                    new_ranks[new_count] = 0;
                    new_index[new_count] = cells.next[cells.next[index] + ich];
                    new_ghost[new_count] = -1;
                    ++new_count;
                }
                continue;
            }

            // Инцидентная ячейка огрубляется. Если узел угловой, тогда роль сохранится,
            // но уже для новой огрубленной ячейки, в остальных случаях запись об
            // инцидентных просто не добавляется.
            if (is_corner_node<dim>(role)) {
                new_roles[new_count] = role;
                new_ranks[new_count] = 0;
                new_index[new_count] = cells.next[index];
                new_ghost[new_count] = -1;
                ++new_count;
                continue;
            }
        }

        for (int j = 0; j < inc_count; ++j) {
            index_t inc = offset + j;
            nodes.incident.role[inc] = new_roles[j];
            nodes.incident.rank[inc] = new_ranks[j];
            nodes.incident.index[inc] = new_index[j];
            nodes.incident.ghost[inc] = new_ghost[j];
        }
    }

    index_t inc_per_node = RawIncident::max_incident_amr(dim);

    for (index_t i = 0; i < split_indices.size(); ++i) {
        index_t ic = split_indices[i];
        int lvl_c = cells.level[ic];

        z_assert(cells.flag[ic] > 0, "Wrong flag");

        index_t in = first_new_node + i * nodes_per_split;
        index_t inc = nodes.incident.offsets[in];

        // Центры дочерних ячеек (4 или 8)
        for (int ich = 0; ich < CpC(dim); ++ich) {
            if constexpr (dim == 2) {
                nodes.incident.role[inc] = SqQuad::iss<0, 0>();
            }
            else {
                nodes.incident.role[inc] = SqCube::iss<0, 0, 0>();
            }
            nodes.incident.rank[inc] = mpi::rank();
            nodes.incident.index[inc] = cells.next[cells.next[ic] + ich];
            nodes.incident.ghost[inc] = -1;

            ++in;
        }

        if constexpr (dim == 2) {
            // +4 (внутренние грани)

            // (i, j) - Индекс дочерней ячейки
            // (I, J) - Роль узла в этой дочерней
            #define set_one(i, j, I, J) { \
                nodes.incident.role[inc] = SqQuad::iss<I, J>(); \
                nodes.incident.rank[inc] = mpi::rank(); \
                nodes.incident.index[inc] = cells.next[cells.next[ic] + Quad::iss<i, j>()]; \
                nodes.incident.ghost[inc] = -1; \
            }

            // Две ячейки у левой грани
            z_assert(inc == nodes.incident.offsets[in + 1], "Bad offset incident #1");
            set_one(-1, -1, 0, +1);
            inc += 1;
            set_one(-1, +1, 0, -1);
            inc += (inc_per_node - 1);

            // Две ячейки у правой грани
            z_assert(inc == nodes.incident.offsets[in + 2], "Bad offset incident #2");
            set_one(+1, -1, 0, +1);
            inc += 1;
            set_one(+1, +1, 0, -1);
            inc += (inc_per_node - 1);

            // Две ячейки у нижней грани
            z_assert(inc == nodes.incident.offsets[in + 3], "Bad offset incident #3");
            set_one(-1, -1, +1, 0);
            inc += 1;
            set_one(+1, -1, -1, 0);
            inc += (inc_per_node - 1);

            // Две ячейки у верхней грани
            z_assert(inc == nodes.incident.offsets[in + 4], "Bad offset incident #4");
            set_one(-1, +1, +1, 0);
            inc += 1;
            set_one(+1, +1, -1, 0);
            inc += (inc_per_node - 1);

            // +8 (внешние, по 2 на каждой грани/на подгранях)
            const SqQuad& orig_cell = cells.verts.mapping<dim>(ic);
            auto orig_children = orig_cell.children();

            for (Side<dim> side: Side<dim>::items()) {
                index_t iface = cells.faces.offsets[ic] + side;
                int symm = cells.faces.adjacent.rotation[iface];

                // Две прилегающих дочерних
                int ich_1 = indexing::children(side)[0];
                int ich_2 = indexing::children(side)[1];

                // индексы вершин грани в ячейках
                int if1 = indexing::amr::sf(side)[0];
                int if2 = indexing::amr::sf(side)[1];

                // координаты новых узлов
                Vector3d v1 = 0.5*(orig_children[ich_1][if1] + orig_children[ich_1][if2]);
                Vector3d v2 = 0.5*(orig_children[ich_2][if1] + orig_children[ich_2][if2]);

                int role_1 = -1, role_2 = -1;
                double eps = 1.0e-5 * (v1 - v2).norm();
                for (int role = 0; role < 9; ++role) {
                    if ((orig_children[ich_1][role] - v1).norm() < eps) {
                        role_1 = role;
                    }
                    if ((orig_children[ich_2][role] - v2).norm() < eps) {
                        role_2 = role;
                    }
                }
                #ifdef SCRUTINY
                if (role_1 < 0 || role_2 < 0) {
                    throw std::runtime_error("Bad find #1944");
                }
                #endif

                if (cells.faces.is_simple(ic, side)) {
                    // сосед того же уровня или ниже
                    z_assert(cells.faces.adjacent.is_local(iface), "Ghosts #123");

                    auto [neibs, jc] = cells.faces.adjacent.get_neib(iface, cells, cells);

                    if (neibs.level[jc] < cells.level[ic]) {
                        // Крупный сосед, должен тоже биться
                        z_assert(neibs.flag[jc] > 0, "Wrong flag #124");

                        // у соседа дочерняя одна прилегает
                        int zch = indexing::adjacent_child(side, symm, cells.z_idx[ic] % CpC(dim));

                        // в этой дочерней две роли
                        int role_1n = -1, role_2n = -1;
                        const SqQuad& child = cells.verts.mapping<dim>(ic);
                        double eps = 1.0e-5 * (v1 - v2).norm();
                        for (int role = 0; role < 9; ++role) {
                            if ((child[role] - v1).norm() < eps) {
                                role_1n = role;
                            }
                            if ((child[role] - v2).norm() < eps) {
                                role_2n = role;
                            }
                        }
                        #ifdef SCRUTINY
                        if (role_1n < 0 || role_2n < 0) {
                            throw std::runtime_error("Bad find #1940");
                        }
                        #endif

                        nodes.incident.role[inc] = role_1n;
                        nodes.incident.rank[inc] = mpi::rank();
                        nodes.incident.index[inc] = cells.next[cells.next[ic] + zch];
                        nodes.incident.ghost[inc] = -1;
                        inc += 1;

                        nodes.incident.role[inc] = role_2n;
                        nodes.incident.rank[inc] = mpi::rank();
                        nodes.incident.index[inc] = cells.next[cells.next[ic] + zch];
                        nodes.incident.ghost[inc] = -1;
                        inc += 1;
                    }


                }
                else {

                }
            }



        }
        else {
            // +12 (центры внутренних граней, по 2 инцидентных)

            // +6 (ребра внутренних, по 4 инцидентных)

            // +24 (внешние подграни)

            // +24 (внешневнутренние ребра, две инцидентных внутри)

            // +24 (внешние ребра, одна инцидентная внутри)

            throw std::runtime_error("Not implemented");
        }
    }
}

} // namespace zephyr::mesh::amr