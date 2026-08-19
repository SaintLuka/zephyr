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

void AmrIncident::reserve(index_t n_nodes, index_t n_values) {
    offsets.reserve(n_nodes + 1);

    role.reserve(n_values);
    rank.reserve(n_values);
    index.reserve(n_values);
    ghost.reserve(n_values);
}

void AmrIncident::shrink_to_fit() {
    offsets.shrink_to_fit();
    role.shrink_to_fit();
    rank.shrink_to_fit();
    index.shrink_to_fit();
    ghost.shrink_to_fit();
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
    coords.clear();
    incident.clear();
}

void AmrNodes::shrink_to_fit() {
    rank.shrink_to_fit();
    next.shrink_to_fit();
    index.shrink_to_fit();
    coords.shrink_to_fit();
    incident.shrink_to_fit();
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
                cells.faces.adjacent.is_alien(iface)) {
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

void AmrNodes::setup_for(AmrCells& cells) {
    // Стираем существующий массив узлов
    clear();

    if (cells.empty()) return;

    // Индексы узлов: index = -13, ghost = -1
    cells.verts.init_unique(-13, -1);

    // TODO: Заменить set, vector на быстрые версии на стеке
    // TODO: Есть только наметки, как это сделать параллельно

    // Последовательная версия работает за один проход по ячейкам
    index_t n_nodes_approx = nodes_estimation(cells.size(), cells.dim());

    coords.reserve(n_nodes_approx);
    rank.reserve(n_nodes_approx);
    next.reserve(n_nodes_approx);
    index.reserve(n_nodes_approx);

    index_t n_inc_approx = cells.dim() < 3 ? 5 * n_nodes_approx : 10 * n_nodes_approx;
    incident.reserve(n_nodes_approx, n_inc_approx);

    int counter = 0;
    for (index_t ic = 0; ic < cells.n_cells(); ++ic) {
        for (role_t loc_iv = 0; loc_iv < cells.verts.max_count(ic); ++loc_iv) {
            index_t iv = cells.verts.offsets[ic] + loc_iv;

            // Нас интересуют актуальные (отмеченные) узлы, которые
            // ещё не получили уникальный индекс.
            if (cells.verts.index[iv] != -13) continue;

            // Ищем все ячейки, которые содержат узел
            auto owners = find_owners(cells, ic, loc_iv);

            // Добавляем узел в массив
            coords.push_back(cells.verts[iv]);
            rank.push_back(0);
            next.push_back(-1);
            index.push_back(counter);

            int n_incident = owners.size();
            if (cells.adaptive()) {
                // Для адаптивных строго фиксируется число инцидентных
                n_incident = cells.dim() == 2 ? AmrIncident::max_amr_count_2D : AmrIncident::max_amr_count_3D;
            }
            incident.offsets.push_back(incident.offsets.back() + n_incident);

            // Отмечаем индекс узла у каждого владельца
            for (auto [ic2, iv2, loc_iv2]: owners) {
                cells.verts.index[iv2] = counter;

                incident.role.push_back(loc_iv2);
                incident.rank.push_back(0);
                incident.index.push_back(ic2);
                incident.ghost.push_back(-1);
            }
            for (int i = owners.size(); i < n_incident; ++i) {
                incident.role.push_back(-1);
                incident.rank.push_back(0);
                incident.index.push_back(-1);
                incident.ghost.push_back(-1);
            }
            ++counter;
        }
    }

    shrink_to_fit();
}

memory_t AmrNodes::memory_usage() const {
    memory_t mem;
    mem.add(next);
    mem.add(rank);
    mem.add(index);
    mem.add(coords);
    return mem;
}

int AmrNodes::check_nodes(const AmrCells& locals) const {
    AmrCells ghosts;
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
    if (coords.size() != n_nodes) {
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
            if (!std::isfinite(coords[in][d])) {
                std::cout << "\tNot finite " << in << " node: " << coords[in].transpose() << "\n";
                return -1;
            }
        }
        if (dim < 3 && coords[in].z() != 0.0) {
            std::cout << "\tNot zero " << in << " node: " << coords[in].transpose() << "\n";
            return -1;
        }
    }

    if (locals.verts.index.size() != locals.verts.size()) {
        std::cout << "\tUnique nodes: bad verts.index.size\n";
        return -1;
    }
    if (locals.verts.ghost.size() != locals.verts.size()) {
        std::cout << "\tUnique nodes: bad verts.ghost.size\n";
        return -1;
    }

    // Проверяем, что в locals.verts верные индексы и координаты
    for (index_t i = 0; i < locals.verts.size(); ++i) {
        Vector3d v1 = locals.verts[i];
        Vector3d v2;
        index_t idx = locals.verts.index[i];
        index_t gst = locals.verts.ghost[i];
        if (gst < 0) {
            // local node
            if (idx < 0 || idx >= n_nodes()) {
                std::cout << "\tLocal vertex index out of range " << idx << " #1\n";
                return -1;
            }
            v2 = coords[idx];
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
            v2 = ghost_nodes.coords[gst];
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

        Vector3d v = coords[in];

        for (index_t i: incident.range(in)) {
            auto role = incident.role[i];
            auto rnk = incident.rank[i];
            auto idx = incident.index[i];
            auto gst = incident.ghost[i];

            if (role < 0) {
                // Фиктивная запись, нет ячейки, только для адаптивных
                if (!locals.adaptive() || idx >= 0 || gst >= 0 || rnk != mpi::rank()) {
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