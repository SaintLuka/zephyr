#include <set>

#include <zephyr/mesh/euler/amr_nodes.h>
#include <zephyr/mesh/euler/amr_cells.h>
#include <zephyr/utils/threads.h>

using zephyr::utils::threads;
using zephyr::geom::Vector3d;
using zephyr::geom::Boundary;

namespace zephyr::mesh {

void AmrIncident::clear() {
    offsets.resize(1);
    offsets[0] = 0;

    role.clear();
    rank.clear();
    index.clear();
    ghost.clear();
}

void AmrIncident::resize(index_t n_values) {
    offsets.resize(n_values + 1, offsets.back());
    role.resize(n_values, -1);
    rank.resize(n_values, -1);
    index.resize(n_values, -1);
    ghost.resize(n_values, -1);
}

void AmrIncident::reserve(index_t n_values) {
    offsets.reserve(n_values + 1);
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
    resize(n_nodes * inc_per_node);
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

inline int nodes_estimation(int n_cells, int dim) {
    z_assert(dim == 2 || dim == 3, "bad dimension");
    if (dim < 3) {
        int nx = int(std::ceil(std::sqrt(n_cells))) + 2;
        return nx * nx;
    } else {
        int nx = int(std::ceil(std::cbrt(n_cells))) + 2;
        return nx * nx * nx;
    }
}

// Владелец узла (любая ячейка, которая содержит узел)
// Основной владелец: ячейка с минимальным индексом.
struct NodeOwner {
    index_t ic; // Индекс ячейки
    index_t iv; // Индекс вершины

    // Считаем различие только по индексу ячейки
    bool operator<(const NodeOwner& other) const { return ic < other.ic; }
};

// Множество ячеек, которые владеют некоторым узлом
struct NodeOwners {
    // Добавить владельца
    void insert(index_t ic, int iv) {
        owners.insert(NodeOwner{ic, iv});
    }

    // Имеется владелец с индексом ic?
    bool contains(index_t ic) const {
        return owners.contains(NodeOwner{ic, -1});
    }

    auto begin() const { return owners.begin(); }

    auto end() const { return owners.end(); }

    std::set<NodeOwner> owners;
};

NodeOwners find_owners(const AmrCells& cells, index_t ic_start, int iv_start) {
    // Интересующая нас вершина
    Vector3d p = cells.verts[iv_start];
    double eps = 1.0e-10 * cells.linear_size(ic_start);

    // Моделирует стек с ячейками в работе
    std::vector<NodeOwner> in_work;

    in_work.emplace_back(NodeOwner{ic_start, iv_start});

    NodeOwners owners;
    while (!in_work.empty()) {
        // Извлекли из стека последнюю ячейку
        auto[ic, iv] = in_work.back();
        in_work.pop_back();

        // Добавили нового владельца
        owners.insert(ic, iv);

        // Проходим по граням, ищем грани, которые содержат искомую вершину.
        // Сосед через такую грань также содержит искомую вершину.
        for (auto iface: cells.faces.range(ic)) {
            if (cells.faces.is_undefined(iface) ||
                cells.faces.is_boundary(iface) ||
                cells.faces.boundary[iface] == Boundary::PERIODIC ||
                cells.faces.adjacent.is_alien(iface)) {
                continue;
            }

            // Проверить, что грань содержит искомую вершину
            bool contain = false;
            for (int j = 0; j < AmrFaces::max_vertices; ++j) {
                int loc_iv = cells.faces.vertices[iface][j];
                if (loc_iv < 0) break;
                if (cells.verts.offsets[ic] + loc_iv == iv) {
                    contain = true;
                    break;
                }
            }

            if (!contain) continue;

            // Грань содержит целевую вершину

            // Индекс соседа через грань
            index_t ic_n = cells.faces.adjacent.index[iface];
            z_assert(ic_n < cells.size(), "Find owners: Out of range");

            // Сосед уже есть в массиве
            if (owners.contains(ic_n)) continue;

            // Ищем интересующую вершину среди вершин соседа
            index_t iv_n = -1;
            for (auto iv2: cells.verts.range(ic_n)) {
                if ((cells.verts[iv2] - p).norm() < eps) {
                    // Проверка на -13??
                    iv_n = iv2;
                    break;
                }
            }

            z_assert(iv_n >= 0, "AmrFaces::setup_for: Impossible error #1");

            // Соседняя ячейка нам подходит, помещаем в стек
            in_work.emplace_back(NodeOwner{ic_n, iv_n});
        }
    }

    return owners;
}

void AmrNodes::setup_for(AmrCells& cells) {
    // Стираем существующий массив узлов
    clear();

    if (cells.empty()) return;

    // Выставляем все индексы на -1
    cells.verts.init_unique();

    // Помечаем актуальные узлы, которые есть на каких-либо гранях, индексом -13.
    threads::parallel_for(
        index_t{0}, cells.n_cells(),
        [&cells](index_t ic) {
            for (auto iface: cells.faces.range(ic)) {
                if (cells.faces.is_undefined(iface)) continue;

                for (int j = 0; j < AmrFaces::max_vertices; ++j) {
                    int loc_iv =  cells.faces.vertices[iface][j];
                    if (loc_iv < 0) break;
                    cells.verts.index[cells.verts.offsets[ic] + loc_iv] = -13;
                }
            }
        });


    // TODO: Заменить set, vector на быстрые версии на стеке
    // TODO: Есть только наметки, как это сделать параллельно

    // Последовательная версия работает за один проход по ячейкам
    /*
     incident:
        std::vector<index_t> offsets = {0};
        std::vector<short> role;
        std::vector<int> rank;
        std::vector<index_t> index;
        std::vector<index_t> alien;
    */

    coords.reserve(nodes_estimation(cells.size(), cells.dim()));

    int counter = 0;
    for (index_t ic = 0; ic < cells.n_cells(); ++ic) {
        for (index_t iv: cells.verts.range(ic)) {
            // Нас интересуют актуальные (отмеченные) узлы, которые
            // ещё не получили уникальный индекс.
            if (cells.verts.index[iv] != -13) continue;

            // Добавляем узел в массив
            coords.push_back(cells.verts[iv]);

            // Ищем все ячейки, которые содержат узел
            auto owners = find_owners(cells, ic, iv);

            // Отмечаем индекс узла у каждого владельца
            for (auto [ic2, iv2]: owners) {
                cells.verts.index[iv2] = counter;
            }
            ++counter;
        }
    }
}

memory_t AmrNodes::memory_usage() const {
    throw std::runtime_error("AmrNodes memory usage not implemented");
}

} // namespace zephyr::mesh