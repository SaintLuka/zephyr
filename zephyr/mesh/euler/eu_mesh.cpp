#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/mesh/primitives/amr_cell.h>
#include <zephyr/mesh/primitives/bfaces.h>
#include <zephyr/mesh/euler/eu_cell.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/utils/json.h>

using zephyr::utils::Json;

namespace zephyr::mesh {

EuMesh::EuMesh(const Json& config, int datasize)
    : m_locals(0, true, datasize), m_aliens(0, true, datasize) {

    if (mpi::master()) {
        // Задать генератор
        auto gen = Generator::create(config);

        initialize(gen->make());
    }

    m_max_level = 0;
    if (config["max_level"]) {
        m_max_level = std::max(0, config["max_level"].as<int>());
        if (m_max_level > 15) {
            std::cerr << "Max level set up to " << m_max_level << ", decreased to 15\n";
            m_max_level = 15;
        }
    }
    if (config["adaptive"]) {
        // Есть ключевое слово adaptive, выставлено на false
        if (!config["adaptive"].as<bool>()) {
            m_max_level = 0;
        }
    }

    // Если есть декомпозиция, то задать задать
    if (!mpi::single()) {
        geom::Box domain = bbox();

        if (config["decomp"]) {
            // Там все свойства декомпозиции автоматически ставятся
            m_decomp = Decomposition::create(domain, config["decomp"]);

            build_aliens();
            redistribute();
        }
        else {
            // По умолчанию что-то такое
            set_decomposition("XY");
        }
    }
}

void EuMesh::initialize(const Grid& grid) {
    m_locals.resize(grid.n_cells());

    for (int i = 0; i < grid.n_cells(); ++i) {
        m_locals[i] = grid.amr_cell(i);
    }

    structured = grid.is_structured();
    m_nx = grid.nx();
    m_ny = grid.ny();
    m_nz = grid.nz();
}

EuCell EuMesh::operator()(int i, int j) {
    i = (i + m_nx) % m_nx;
    j = (j + m_ny) % m_ny;
    return operator[](m_ny * i + j);
}

EuCell EuMesh::operator()(int i, int j, int k) {
    i = (i + m_nx) % m_nx;
    j = (j + m_ny) % m_ny;
    k = (k + m_nz) % m_nz;
    return operator[](m_nz * (m_ny * i + j) + k);
}

geom::Box EuMesh::bbox() {
    geom::Box box1 = Box::Empty(3);
    for (auto& cell: *this) {
        for (auto& face: cell.faces()) {
            for (int i = 0; i < face.size(); ++i) {
                box1.capture(face.vs(i));
            }
        }
    }

    geom::Box box2;

#ifdef ZEPHYR_MPI
    if (!mpi::single()) {
        // Покомпонентный минимум/максимум
        MPI_Allreduce(box1.vmin.data(), box2.vmin.data(), 3, MPI_DOUBLE, MPI_MIN, mpi::comm());
        MPI_Allreduce(box1.vmax.data(), box2.vmax.data(), 3, MPI_DOUBLE, MPI_MAX, mpi::comm());
    }
#endif

    return box2;
}

bool EuMesh::has_nodes() const {
    return !m_nodes.empty();
}

void EuMesh::break_nodes() {
    m_nodes.clear();
}

inline int nodes_estimation(int n_cells, int dim) {
    assert(dim == 2 || dim == 3);
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
    int ic;  // Индекс ячейки
    int iv;  // Индекс вершины в ячейке

    // Можно создать без указания iv, типа просто индекс ячейки
    // используется для поиска в множестве std::set<NodeOwner>;
    NodeOwner(int _ic, int _iv)
            : ic(_ic), iv(_iv) {}

    // Считаем различие только по индексу ячейки
    inline bool operator<(const NodeOwner& other) const {
        return ic < other.ic;
    }

    static const int base = 64;

    // Смешаный индекс вершины (владелец + индекс внутри)
    inline int node_index() const {
        return base * ic + iv;
    }

    // Восстановление из смешаного индекса
    inline static NodeOwner from_index(int idx) {
        return NodeOwner(idx / base, idx % base);
    }
};

// Простая структура, хранит множество ячеек, которые владеют
// некоторым узлом
struct NodeOwners {

    // Добавить владельца
    inline void insert(int ic, int iv) {
        owners.insert(NodeOwner(ic, iv));
    }

    // Имеется владелец с индексом ic?
    inline bool contain(int ic) {
        return owners.count(NodeOwner(ic, -1)) > 0;
    }

    inline NodeOwner main() const {
        return *std::min_element(owners.begin(), owners.end());
    }

    inline auto begin() {
        return owners.begin();
    }

    inline auto end() {
        return owners.end();
    }

    std::set<NodeOwner> owners;
};

NodeOwners find_owners(AmrStorage& cells, AmrCell& base_cell, int base_iv) {
    NodeOwners owners;

    // Интересующая нас вершина
    Vector3d base_v = base_cell.vertices[base_iv];
    double eps = 1.0e-10 * base_cell.size;

    // Моделирует стек с ячейками в работе
    std::vector<NodeOwner> in_work;

    in_work.emplace_back(base_cell.index, base_iv);

    while (!in_work.empty()) {
        // Извлекли из стека последнюю ячейку
        auto[ic_c, iv_c] = in_work.back();
        in_work.pop_back();

        owners.insert(ic_c, iv_c);

        AmrCell &cell = cells[ic_c];

        // Проходим по граням, ищем грани, которые содержат искомую вершину
        // и при этом указывают на ячейки, которых нет в множестве done
        for (const BFace &face: cell.faces) {
            if (face.is_undefined() ||
                face.is_boundary()) {
                continue;
            }

            // Грань содержит целевую вершину
            if (face.contain(iv_c)) {
                // Сосед через грань
                AmrCell &neib = cells[face.adjacent.index];
                int ic_n = neib.index;

                // Сосед уже есть в массиве
                if (owners.contain(ic_n)) {
                    continue;
                }

                // Ищем интересующую вершину среди вершин соседа
                int iv_n = neib.vertices.find(base_v, eps);

                assert(iv_n >= 0 && "EuMesh::collect_nodes: Impossible error #1");

                // Соседняя ячейка нам подходит, помещаем в стек
                in_work.emplace_back(ic_n, iv_n);
            }
        }
    }

    return owners;
}

// Для узлов ячейки, которые помечены индексом -13 выставляет их основного
// владельца. Кроме этого возвращает количество узлов, для которых ячейка
// является основным владельцем.
int count_owned(AmrStorage::Item& cell, AmrStorage& cells) {
    int owned = 0;
    for (int iv = 0; iv < BNodes::max_count; ++iv) {
        // Интересуют актуальные (помеченые) узлы,
        // которые ещё не получили свой номер.
        if (cell.nodes[iv] != -13) {
            continue;
        }

        // Ищем все ячейки, которые содержат узел
        auto owners = find_owners(cells, cell, iv);

        // Основной владелец узла
        NodeOwner owner = owners.main();

        // Устанавливаем смешаный индекс (кодирует основного
        // владельца) для каждого узла
        cell.nodes[iv] = owner.node_index();

        if (owner.ic == cell.index) {
            ++owned;
        }
    }
    return owned;
};

void setup_nodes(AmrStorage& cells, std::vector<Vector3d>& nodes) {
    if (cells.empty()) {
        return;
    }

    // Стираем существующий массив узлов
    nodes.clear();

    // Помечаем актуальные узлы (которые есть на каких-либо гранях),
    // индексом -13.
    threads::for_each(cells.begin(), cells.end(),
            [](AmrStorage::Item& cell) {
                cell.mark_actual_nodes(-13);
            });

    // TODO: Заменить set, vector на быстрые версии на стеке

    // TODO: Есть только наметки, как это сделать параллельно

    if (true || !threads::active()) {
        /// Последовательная версия работает за один проход по ячейкам
        nodes.reserve(nodes_estimation(cells.size(), cells[0].dim));

        int counter = 0;
        for (auto &cell: cells) {
            for (int iv = 0; iv < BNodes::max_count; ++iv) {
                // Интересуют актуальные (помеченые) узлы,
                // которые ещё не получили свой номер.
                if (cell.nodes[iv] != -13) {
                    continue;
                }

                // Ищем все ячейки, которые содержат узел
                auto owners = find_owners(cells, cell, iv);

                // Отмечаем индекс узла у каждого владельца
                for (auto &owner: owners) {
                    cells[owner.ic].nodes[owner.iv] = counter;
                }

                // Добавляем узел в массив
                nodes.push_back(cell.vertices[iv]);
                ++counter;
            }
        }
    }
    else {
        // Параллельная версия отличается от последовательной и
        // содержит две стадии.

        std::vector<int> n_owned = threads::partial_sum(
                cells.begin(), cells.end(), 0,
                count_owned, std::ref(cells));

        // Смещение индексации в каждом треде
        std::vector<int> offset;
        offset.reserve(n_owned.size());
        offset.push_back(0);
        for (size_t i = 1; i < n_owned.size(); ++i) {
            offset[i] = offset[i - 1] + n_owned[i];
        }

        int counter = 0;
        for (auto &cell: cells) {

            // Идем по всем вершинам с выставленным номером
            for (int iv = 0; iv < BNodes::max_count; ++iv) {
                if (cell.nodes[iv] < 0) {
                    continue;
                }

                NodeOwner owner = NodeOwner::from_index(cell.nodes[iv]);

                if (owner.ic == cell.index) {
                    // Вершина принадлежит ячейке
                    nodes.push_back(cell.vertices[iv]);

                    cell.nodes[iv] = counter;
                    ++counter;
                }
                else {
                    assert(owner.ic < cell.index && "setup unique, error #1534");
                    cell.nodes[iv] = cells[owner.ic].nodes[owner.iv];
                }


            }
        }
    }
}

void EuMesh::collect_nodes() {
    if (has_nodes()) {
        return;
    }
    setup_nodes(m_locals, m_nodes);
}

} // namespace zephyr::mesh