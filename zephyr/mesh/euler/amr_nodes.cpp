#include <set>

#include <zephyr/utils/mpi.h>
#include <zephyr/mesh/euler/amr_nodes.h>
#include <zephyr/mesh/euler/amr_cells.h>
#include <zephyr/utils/threads.h>

using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::geom::Vector3d;
using zephyr::geom::Boundary;

namespace zephyr::mesh {

void AmrIncident::clear() {
    offsets = {0};

    role.clear();
    rank.clear();
    index.clear();
    ghost.clear();
}

void AmrIncident::resize(index_t n_nodes, index_t n_values) {
    offsets.resize(n_nodes + 1, offsets.back());

    role.resize(n_values, -1);
    rank.resize(n_values, -1);
    index.resize(n_values, -1);
    ghost.resize(n_values, -1);
}

void AmrIncident::resize_amr(index_t n_nodes, int dim) {
    z_assert(dim == 2 || dim == 3, "AmrIncident::resize: bad dimension");

    index_t prev_size = role.size();
    int inc_per_node = dim < 3 ? max_amr_count_2D : max_amr_count_3D;
    resize(n_nodes, n_nodes * inc_per_node);
    for (index_t i = prev_size; i < role.size(); ++i) {
        offsets[i + 1] = offsets[i] + inc_per_node;
    }
}

void AmrIncident::reserve(index_t n_nodes, index_t n_values) {
    offsets.reserve(n_nodes + 1);

    role.reserve(n_values);
    rank.reserve(n_values);
    index.reserve(n_values);
    ghost.reserve(n_values);
}

void AmrIncident::reserve_amr(index_t n_nodes, int dim) {
    z_assert(dim == 2 || dim == 3, "AmrIncident::reserve: bad dimension");

    int inc_per_node = max_amr_count(dim);
    reserve(n_nodes, inc_per_node * n_nodes);
}

void AmrIncident::shrink_to_fit() {
    offsets.shrink_to_fit();
    role.shrink_to_fit();
    rank.shrink_to_fit();
    index.shrink_to_fit();
    ghost.shrink_to_fit();
}

int AmrIncident::count(index_t inode) const {
    index_t i = offsets[inode];
    while (i < offsets[inode + 1] && role[i] >= 0) {
        ++i;
    }
    return i - offsets[inode];
}

memory_t AmrIncident::memory_usage() const {
    memory_t mem;
    mem.add(role);
    mem.add(rank);
    mem.add(index);
    mem.add(ghost);
    return mem;
}

void AmrNodes::clear() {
    rank.clear();
    next.clear();
    index.clear();
    coord.clear();
    incident.clear();
}

void AmrNodes::resize(index_t n_nodes, index_t n_incident) {
    data.resize(n_nodes);
    rank.resize(n_nodes);
    next.resize(n_nodes);
    index.resize(n_nodes);
    coord.resize(n_nodes);
    incident.resize(n_nodes, n_incident);
}

void AmrNodes::resize_amr(index_t n_nodes, int dim) {
    data.resize(n_nodes);
    rank.resize(n_nodes);
    next.resize(n_nodes);
    index.resize(n_nodes);
    coord.resize(n_nodes);
    incident.resize_amr(n_nodes, dim);
}

void AmrNodes::reserve(index_t n_nodes, index_t n_incident) {
    data.reserve(n_nodes);
    rank.reserve(n_nodes);
    next.reserve(n_nodes);
    index.reserve(n_nodes);
    coord.reserve(n_nodes);
    incident.reserve(n_nodes, n_incident);
}

void AmrNodes::reserve_amr(index_t n_nodes, int dim) {
    data.reserve(n_nodes);
    rank.reserve(n_nodes);
    next.reserve(n_nodes);
    index.reserve(n_nodes);
    coord.reserve(n_nodes);
    incident.reserve_amr(n_nodes, dim);
}

void AmrNodes::shrink_to_fit() {
    data.shrink_to_fit();
    rank.shrink_to_fit();
    next.shrink_to_fit();
    index.shrink_to_fit();
    coord.shrink_to_fit();
    incident.shrink_to_fit();
}

void AmrNodes::copy_data(index_t from, index_t to) {
    copy_data(from, this, to);
}

void AmrNodes::copy_data(index_t from, AmrNodes* dst, index_t to) const {
    data.copy_data(from, &dst->data, to);
}

void AmrNodes::copy_geom(index_t in, AmrNodes& nodes, index_t jn, index_t inc_offset) const {
    nodes.rank [jn] = rank [in];
    nodes.next [jn] = next [in];
    nodes.index[jn] = index[in];
    nodes.coord[jn] = coord[in];

    nodes.incident.offsets[jn] = inc_offset;
    nodes.incident.offsets[jn + 1] = inc_offset + incident.max_count(in);

    for (index_t i = 0; i < incident.max_count(in); ++i) {
        index_t I = incident.offsets[in] + i;
        index_t J = nodes.incident.offsets[jn] + i;

        nodes.incident.role [J] = incident.role [I];
        nodes.incident.rank [J] = incident.rank [I];
        nodes.incident.index[J] = incident.index[I];
        nodes.incident.ghost[J] = incident.ghost[I];
    }
}

inline index_t nodes_estimation(int n_cells, int dim) {
    z_assert(dim == 2 || dim == 3, "bad dimension");
    if (dim < 3) {
        int nx = index_t(std::ceil(std::sqrt(n_cells))) + 2;
        return nx * nx;
    } else {
        int nx = index_t(std::ceil(std::cbrt(n_cells))) + 2;
        return nx * nx * nx;
    }
}

// Владелец узла (любая ячейка, которая содержит узел)
// Основной владелец: ячейка с минимальным индексом.
struct NodeOwner {
    index_t ic; // Индекс ячейки
    index_t iv; // Индекс вершины (в массиве дублирующихся verts)
    role_t  role; // Локальный индекс вершины в ячейке

    // Считаем различие только по индексу ячейки
    bool operator<(const NodeOwner& other) const { return ic < other.ic; }
};

// Множество ячеек, которые владеют некоторым узлом
struct NodeOwners {
    // Добавить владельца
    void insert(index_t ic, index_t iv, role_t role) {
        owners.insert(NodeOwner{.ic=ic, .iv=iv, .role=role});
    }

    // Имеется владелец с индексом ic?
    bool contains(index_t ic) const {
        return owners.contains(NodeOwner{.ic=ic, .iv=-1, .role=0});
    }

    auto begin() const { return owners.begin(); }

    auto end() const { return owners.end(); }

    auto size() const { return owners.size(); }

    std::set<NodeOwner> owners;
};

NodeOwners find_owners(const AmrCells& cells, index_t ic0, role_t loc_iv0) {
    // Интересующая нас вершина
    index_t iv0 = cells.verts.offsets[ic0] + loc_iv0;
    Vector3d p = cells.verts[iv0];
    double eps = 1.0e-10 * cells.linear_size(ic0);

    // Моделирует стек с ячейками в работе
    std::vector<NodeOwner> in_work;

    in_work.emplace_back(NodeOwner{.ic=ic0, .iv=iv0, .role=loc_iv0});

    NodeOwners owners;
    while (!in_work.empty()) {
        // Извлекли из стека последнюю ячейку
        auto[ic, iv, loc_iv] = in_work.back();
        in_work.pop_back();

        // Добавили нового владельца
        owners.insert(ic, iv, loc_iv);

        // Ищем грани, которые содержат искомую вершину, сосед
        // через такую грань вероятно также содержит эту вершину.
        for (auto iface: cells.faces.range(ic)) {
            if (cells.faces.is_undefined(iface) ||
                cells.faces.is_boundary(iface) ||
                cells.faces.boundary[iface] == Boundary::PERIODIC ||
                cells.faces.adjacent.is_ghost(iface)) {
                continue;
            }

            // Проверить, что грань содержит целевую вершину
            auto face_n = cells.faces.normal[iface];
            auto face_c = cells.faces.center[iface];

            // Вершина не в плоскости грани
            if (std::abs((p - face_c).dot(face_n)) > eps) continue;

            // Индекс соседа через грань
            index_t ic_n = cells.faces.adjacent.index[iface];
            z_assert(ic_n < cells.size(), "Find owners: Out of range");

            // Сосед уже есть в массиве
            if (owners.contains(ic_n)) continue;

            // Грань содержит целевую вершину

            // Ищем интересующую вершину среди вершин соседа
            index_t iv_n = -1;
            role_t loc_iv_n = 0;
            for (; loc_iv_n < cells.verts.max_count(ic_n); ++loc_iv_n) {
                index_t iv2 = cells.verts.offsets[ic_n] + loc_iv_n;
                if ((cells.verts[iv2] - p).norm() < eps) {
                    iv_n = iv2;
                    break;
                }
            }

            // Сосед может не содержать искомую вершину
            if (loc_iv_n >= cells.verts.max_count(ic_n)) continue;

            z_assert(iv_n >= 0, "AmrFaces::setup_for: Impossible error");

            // Соседняя ячейка нам подходит, помещаем в стек
            in_work.emplace_back(NodeOwner{.ic=ic_n, .iv=iv_n, .role=loc_iv_n});
        }
    }

    return owners;
}

template <bool complete>
AmrNodes::Incomplete AmrNodes::generate(const AmrCells& cells) {
    AmrNodes nodes;
    AmrVerts verts;

    // Пустые массивы
    if (cells.empty()) {
        return {verts, nodes};
    }

    // Индексы узлов: index = -13
    verts.index.resize(cells.n_verts());
    std::ranges::fill(verts.index, -13);

    // TODO: Заменить set, vector на быстрые версии на стеке
    // TODO: Есть только наметки, как это сделать параллельно

    // Последовательная версия работает за один проход по ячейкам
    index_t n_nodes_approx = nodes_estimation(cells.size(), cells.dim());

    nodes.coord.reserve(n_nodes_approx);
    if constexpr (complete) {
        nodes.data.reserve(n_nodes_approx);
        nodes.rank.reserve(n_nodes_approx);
        nodes.next.reserve(n_nodes_approx);
        nodes.index.reserve(n_nodes_approx);

        index_t n_inc_approx = cells.dim() < 3 ? 5 * n_nodes_approx : 10 * n_nodes_approx;
        nodes.incident.reserve(n_nodes_approx, n_inc_approx);
    }

    int counter = 0;
    for (index_t ic = 0; ic < cells.n_cells(); ++ic) {
        for (role_t loc_iv = 0; loc_iv < cells.verts.max_count(ic); ++loc_iv) {
            index_t iv = cells.verts.offsets[ic] + loc_iv;

            // Нас интересуют актуальные (отмеченные) узлы, которые
            // ещё не получили уникальный индекс.
            if (verts.index[iv] != -13) continue;

            // Ищем все ячейки, которые содержат узел
            auto owners = find_owners(cells, ic, loc_iv);

            // Добавляем узел в массив
            nodes.coord.push_back(cells.verts[iv]);

            int n_incident = owners.size();
            if (cells.adaptive()) {
                // Для адаптивных строго фиксируется число инцидентных
                n_incident = cells.dim() == 2 ? AmrIncident::max_amr_count_2D : AmrIncident::max_amr_count_3D;
            }

            if constexpr (complete) {
                nodes.data.resize(nodes.data.size() + 1);
                nodes.rank.push_back(0);
                nodes.next.push_back(-1);
                nodes.index.push_back(counter);

                nodes.incident.offsets.push_back(nodes.incident.offsets.back() + n_incident);
            }

            // Отмечаем индекс узла у каждого владельца
            for (auto [ic2, iv2, loc_iv2]: owners) {
                verts.index[iv2] = counter;

                if constexpr (complete) {
                    nodes.incident.role.push_back(loc_iv2);
                    nodes.incident.rank.push_back(0);
                    nodes.incident.index.push_back(ic2);
                    nodes.incident.ghost.push_back(-1);
                }
            }
            if constexpr (complete) {
                for (int i = owners.size(); i < n_incident; ++i) {
                    nodes.incident.role.push_back(-1);
                    nodes.incident.rank.push_back(0);
                    nodes.incident.index.push_back(-1);
                    nodes.incident.ghost.push_back(-1);
                }
            }
            ++counter;
        }
    }
    nodes.shrink_to_fit();
    return {verts, nodes};
}

template AmrNodes::Incomplete AmrNodes::generate<true >(const AmrCells& cells);
template AmrNodes::Incomplete AmrNodes::generate<false>(const AmrCells& cells);

void AmrNodes::setup_for(AmrCells& cells) {
    auto [verts, nodes] = AmrNodes::generate<true>(cells);

    z_assert(cells.has_nodes(), "No nodes");
    z_assert(cells.n_verts() == verts.index.size(), "bad sizes");
    cells.verts.index = std::move(verts.index);

    *this = std::move(nodes);
}

memory_t AmrNodes::memory_usage() const {
    memory_t mem;
    mem.add(next);
    mem.add(rank);
    mem.add(index);
    mem.add(coord);
    return mem;
}

int AmrNodes::check_nodes(const AmrCells& locals) const {
    AmrCells ghosts = locals.same();
    AmrNodes ghost_nodes;
    return check_nodes(locals, ghosts, ghost_nodes);
}

int AmrNodes::check_sizes() const {
    int n_nodes = size();
    if (n_nodes < 1) {
        std::cout << "\tHas no unique nodes\n";
        return -1;
    }
    if (rank.size() != n_nodes) {
        std::cout << "\tUnique nodes: bad nodes.rank.size\n";
        return -1;
    }
    if (next.size() != n_nodes) {
        std::cout << "\tUnique nodes: bad nodes.next.size\n";
        return -1;
    }
    if (index.size() != n_nodes) {
        std::cout << "\tUnique nodes: bad nodes.index.size\n";
        return -1;
    }
    if (coord.size() != n_nodes) {
        std::cout << "\tUnique nodes: bad nodes.coords.size\n";
        return -1;
    }
    if (incident.offsets.size() != n_nodes + 1) {
        std::cout << "\tUnique nodes: bad nodes.incident.offsets.size\n";
        return -1;
    }
    return 0;
}

int AmrNodes::check_nodes(const AmrCells& locals, const AmrCells& ghosts, const AmrNodes& ghost_nodes) const {
    int res = check_sizes();
    if (res < 0) return res;

    // Базовые проверки
    int dim = locals.dim();
    for (index_t in = 0; in < n_nodes(); ++in) {
        if (index[in] < 0 || index[in] != in) {
            std::cout << "\tWrong node index\n";
            return -1;
        }
        if (rank[in] < 0 || rank[in] != mpi::rank()) {
            std::cout << "\tWrong node rank\n";
            return -1;
        }
        for (int d = 0; d < dim; ++d) {
            if (!std::isfinite(coord[in][d])) {
                std::cout << "\tNot finite " << in << " node: " << coord[in].transpose() << "\n";
                return -1;
            }
        }
        if (dim < 3 && coord[in].z() != 0.0) {
            std::cout << "\tNot zero " << in << " node: " << coord[in].transpose() << "\n";
            return -1;
        }
    }

    if (locals.verts.rank.size() != locals.verts.n_verts()) {
        std::cout << "\tUnique nodes: bad verts.rank.size\n";
        return -1;
    }
    if (locals.verts.index.size() != locals.verts.n_verts()) {
        std::cout << "\tUnique nodes: bad verts.index.size\n";
        return -1;
    }
    if (locals.verts.ghost.size() != locals.verts.n_verts()) {
        std::cout << "\tUnique nodes: bad verts.ghost.size\n";
        return -1;
    }

    // Проверяем, что в locals.verts верные индексы и координаты
    for (index_t i = 0; i < locals.verts.n_verts(); ++i) {
        Vector3d v1 = locals.verts[i];
        Vector3d v2;
        index_t rnk = locals.verts.rank[i];
        index_t idx = locals.verts.index[i];
        index_t gst = locals.verts.ghost[i];
        if (gst < 0) {
            // local node
            if (idx < 0 || idx >= n_nodes()) {
                std::cout << "\tLocal vertex index out of range " << idx << " #1\n";
                return -1;
            }
            if (rnk != mpi::rank() || rnk != rank[idx]) {
                std::cout << "\tLocal vertex bad rank " << rnk << "\n";
                return -1;
            }
            v2 = coord[idx];
        }
        else {
            // ghost node
#ifndef ZEPHYR_MPI
            std::cout << "\tGhost index >= 0 but mpi is not enabled\n";
            return -1;
#else
            if (idx < 0) {
                std::cout << "\tLocal vertex index out of range " << idx << " #2\n";
                return -1;
            }
            if (mpi::single()) {
                std::cout << "\tGhost index >= 0 for single mpi run\n";
                return -1;
            }
            if (gst < 0 || gst >= ghost_nodes.size()) {
                std::cout << "\tGhost vertex index out of range " << gst << "\n";
                return -1;
            }
            if (rnk == mpi::rank() || rnk != ghost_nodes.rank[gst]) {
                std::cout << "\tRemote vertex bad rank " << rnk << "\n";
                std::cout << "\t\tRank: " << mpi::rank() << "; iv: " << i << "\n";
                return -1;
            }
            v2 = ghost_nodes.coord[gst];
#endif
        }
        if (v1 != v2) {
            std::cout << "\tDifferent coords " << v1.transpose() << "; " << v2.transpose() << "; diff: " << (v1 - v2).norm() << "\n";
            return -1;
        }
    }

    // Проверяем смежность
    for (index_t in = 0; in < n_nodes(); ++in) {
        if (locals.adaptive()) {
            if ((locals.dim() == 2 && incident.max_count(in) != AmrIncident::max_amr_count_2D) ||
                (locals.dim() == 3 && incident.max_count(in) != AmrIncident::max_amr_count_3D)) {
                std::cout << "Wrong max count of incident cells: " << incident.max_count(in) << "\n";
                return -1;
            }
        }

        Vector3d v = coord[in];

        for (index_t i: incident.range(in)) {
            auto role = incident.role[i];
            auto rnk = incident.rank[i];
            auto idx = incident.index[i];
            auto gst = incident.ghost[i];

            if (role < 0) {
                // Фиктивная запись, нет ячейки, только для адаптивных
                if (!locals.adaptive() || idx >= 0 || gst >= 0) { // || rnk != mpi::rank()) { пока убрал строгую проверку
                    std::cout << "\tBad fake incident cell: " << idx << ", " << gst << ", " << rnk << "\n";
                    return -1;
                }
            }
            else if (gst < 0) {
                // local incident cell
                if (idx < 0 || idx >= locals.size()) {
                    std::cout << "\tLocal incident cell index out of range " << idx << " #1\n";
                    return -1;
                }
                if (locals.rank[idx] != rnk) {
                    std::cout << "\tBad incident cell rank\n";
                    return -1;
                }
                if (role < 0 || role >= locals.verts.count(idx)) {
                    std::cout << "\tBad incident cell role #1: " << int(role) << "/" << locals.verts.count(idx) << "\n";
                    return -1;
                }
                if (locals.vertex(idx, role) != v) {
                    std::cout << "\tBad incident cell role #2\n";
                    return -1;
                }
            }
            else {
                // ghost node
#ifndef ZEPHYR_MPI
                std::cout << "\tGhost incident cell index >= 0 but mpi is not enabled\n";
                return -1;
#else
                if (idx < 0) {
                    std::cout << "\tLocal incident cell index out of range " << idx << " #2\n";
                    return -1;
                }
                if (mpi::single()) {
                    std::cout << "\tGhost incident cell index >= 0 for single mpi run\n";
                    return -1;
                }
                if (gst < 0 || gst >= ghosts.size()) {
                    std::cout << "\tGhost incident cell index out of range " << gst << "\n";
                    return -1;
                }
                if (ghosts.rank[gst] != rnk) {
                    std::cout << "\tBad incident cell rank (ghost)\n";
                    return -1;
                }
                if (role < 0 || role >= ghosts.verts.count(gst)) {
                    std::cout << "\tBad incident cell role (ghost) #1\n";
                    return -1;
                }
                if (ghosts.vertex(gst, role) != v) {
                    std::cout << "\tBad incident cell role (ghost) #2\n";
                    return -1;
                }
#endif
            }
        }
    }
    return 0;
}

} // namespace zephyr::mesh