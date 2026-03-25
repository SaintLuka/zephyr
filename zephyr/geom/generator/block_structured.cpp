#include <array>
#include <filesystem>
#include <algorithm>
#include <fstream>
#include <format>
#include <map>
#include <list>
#include <ranges>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/side.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>
#include <zephyr/geom/generator/curve/plane.h>
#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/block_structured.h>
#include <zephyr/utils/threads.h>

namespace zephyr::geom::generator {

inline bool good_number(double num) {
    return !std::isinf(num) && !std::isnan(num) && num > 0.0;
}

inline bool bad_number(double num) {
    return std::isinf(num) || std::isnan(num) || num <= 0.0;
}

constexpr int node_idx(Side side, int idx) {
    z_assert(0 <= static_cast<int>(side) && static_cast<int>(side) < 4,
        "node_idx: invalid side");
    switch (side) {
        case Side::L: return std::array{0, 2}[idx];
        case Side::R: return std::array{1, 3}[idx];
        case Side::B: return std::array{0, 1}[idx];
        case Side::T:
        default:      return std::array{2, 3}[idx];
    }
}

// Перпендикулярная ось
constexpr Axis orthogonal(Axis axis) {
    return axis == Axis::X ? Axis::Y : Axis::X;
}

// Противоположная сторона
constexpr Side opposite(Side side) {
    z_assert(0 <= static_cast<int>(side) && static_cast<int>(side) < 4,
        "opposite: invalid side");
    switch (side) {
        case Side::L: return Side::R;
        case Side::R: return Side::L;
        case Side::B: return Side::T;
        case Side::T: ;
        default:      return Side::B;
    }
}

// Стороны, которые лежат вдоль выбранной оси
constexpr std::array<Side, 2> sides_by_axis(Axis axis) {
    return axis == Axis::X? std::array{Side::B, Side::T} : std::array{Side::L, Side::R};
}

// Пара сторон, связанных с вершиной (обход против часовой внутри блока)
constexpr std::tuple<Side, Side> node_incident_sides(int v_idx) {
    z_assert(0 <= v_idx && v_idx < 4, "Invalid base node index");
    switch (v_idx) {
        case 0: return {Side::B, Side::L};
        case 1: return {Side::R, Side::B};
        case 2: return {Side::L, Side::T};
        case 3:
        default: return {Side::T, Side::R};
    }
}

// ---------- Самые простые конструкторы, get/set методы ------------- -------------------------------------------------

BlockStructured::BlockStructured()
    : Generator("block-structured") {
}

int BlockStructured::n_blocks() const {
    return static_cast<int>(m_blocks.size());
}

void BlockStructured::disable_verbosity() {
    m_verbosity = 0;
}

void BlockStructured::set_verbosity(int verbose) {
    m_verbosity = std::max(0, verbose);
}

void BlockStructured::set_iters_count(int iters_count) {
    m_iters_count = std::max(0, iters_count);
}

// ---------- Редактирование ----------------------------------------- -------------------------------------------------

Block &BlockStructured::operator[](int idx) const {
    if (m_stage != EDITABLE) {
        throw std::runtime_error("BlockStructured::operator[] is only available at the editing stage");
    }
    if (idx >= n_blocks() || idx < -n_blocks()) {
        throw std::out_of_range("BlockStructured::operator[]: index > blocks.size()");
    }
    if (idx < 0) {
        idx += n_blocks();
    }
    if (!m_blocks[idx]) {
        throw std::runtime_error("BlockStructured::operator[]: block does not exist");
    }
    return *m_blocks[idx];
}

inline void skip_warning(const std::string& func_name) {
    std::cerr << "Skip " << func_name << ". It is only available at the editing stage\n";
}

void BlockStructured::operator+=(const std::array<BaseNode::Ptr, 4>& nodes) {
    if (m_stage != EDITABLE) {
        skip_warning("BlockStructured::operator+=()");
        return;
    }
    Block::Ptr block = Block::create(nodes);
    block->set_index(n_blocks());
    m_blocks.emplace_back(block);

    // Добавить рёбра
    for (Side side: sides_2D) {
        auto v1 = block->base_node(side, 0);
        auto v2 = block->base_node(side, 1);
        BaseEdge edge{v1, v2};

        if (m_edges.contains(edge)) {
            m_edges[edge].add(block, side);
        }
        else {
            m_edges[edge] = BlockPair{.b1=block, .side1=side};
        }
    }
}

Curve::Ptr BlockStructured::get_boundary(BaseNode::Ref v1, BaseNode::Ref v2) {
    BaseEdge edge{v1, v2};
    if (m_edges.contains(edge)) {
        if (m_edges[edge].empty()) {
            std::cerr << "BlockStructured::get_boundary warning: has no blocks on edge\n";
            return nullptr;
        }
        else if (m_edges[edge].boundary()) {
            return m_edges[edge].b1.lock()->boundary(m_edges[edge].side1);
        }
        else {
            std::cerr << "BlockStructured::get_boundary warning: attempt to get boundary for inner edge\n";
            return nullptr;
        }
    }
    else {
        throw std::runtime_error("BlockStructured::get_boundary: has no edge for pair of base nodes");
    }
}

void BlockStructured::set_boundary(BaseNode::Ref v1, BaseNode::Ref v2, Curve::Ref curve) {
    if (m_stage != EDITABLE) {
        skip_warning("BlockStructured::set_boundary");
        return;
    }
    BaseEdge edge{v1, v2};
    if (m_edges.contains(edge)) {
        if (m_edges[edge].empty()) {
            std::cerr << "BlockStructured::set_boundary warning: has no blocks on edge\n";
        }
        else if (m_edges[edge].boundary()) {
            m_edges[edge].b1.lock()->set_boundary(m_edges[edge].side1, curve);
        }
        else {
            std::cerr << "BlockStructured::set_boundary warning: attempt to set boundary for inner edge\n";
        }
    }
    else {
        throw std::runtime_error("BlockStructured::set_boundary: has no edge for pair of base nodes");
    }
}

void BlockStructured::set_boundary(BaseNode::Ref v1, BaseNode::Ref v2) {
    set_boundary(v1, v2, Boundary::UNDEFINED);
}

void BlockStructured::set_boundary(BaseNode::Ref v1, BaseNode::Ref v2, Boundary boundary) {
    if (m_stage != EDITABLE) {
        skip_warning("BlockStructured::set_boundary");
        return;
    }
    BaseEdge edge{v1, v2};
    if (m_edges.contains(edge)) {
        if (m_edges[edge].empty()) {
            std::cerr << "BlockStructured::set_boundary warning: has no blocks on edge\n";
            return;
        }
        else if (m_edges[edge].boundary()) {
            auto block = m_edges[edge].b1.lock();
            auto side = m_edges[edge].side1;
            if (block->adjacent_block(side)) {
                std::cerr << "BlockStructured::set_boundary warning: can't set boundary to inner edge\n";
                return;
            }
            if (block->boundary(side)) {
                if (boundary != Boundary::UNDEFINED) {
                    block->boundary(side)->set_boundary(boundary);
                }
                else {
                    std::cerr << "BlockStructured::set_boundary warning: boundary curve is already defined\n";
                }
            }
            else {
                Curve::Ptr plane = Plane::create(v1, v2);
                if (boundary != Boundary::UNDEFINED) {
                    plane->set_boundary(boundary);
                }
                block->set_boundary(side, plane);
            }
        }
        else {
            std::cerr << "BlockStructured::set_boundary warning: attempt to set boundary for inner edge\n";
            return;
        }
    }
    else {
        throw std::runtime_error("BlockStructured::set_boundary: has no edge for pair of base nodes");
    }
}

void BlockStructured::set_boundary(const std::vector<BaseNode::Ptr>& nodes, Curve::Ref curve) {
    if (m_stage != EDITABLE) {
        skip_warning("BlockStructured::set_boundary");
        return;
    }
    if (nodes.size() < 2) {
        throw std::invalid_argument("BlockStructured::set_boundary: at least two nodes required");
    }
    for (int i = 1; i < nodes.size(); ++i) {
        set_boundary(nodes[i - 1], nodes[i], curve);
    }
}

// Точки лежат на одной прямой или нет?
bool collinear(const std::vector<BaseNode::Ptr>& nodes, double eps = 1e-4) {
    z_assert(nodes.size() >= 2, "collinear: base size");

    if (nodes.size() == 2) return true;

    // Берём первую и последнюю точку для построения прямой
    Vector3d v1 = nodes.front()->pos();
    Vector3d v2 = nodes.back()->pos();
    Vector3d c = 0.5 * (v1 + v2);
    Vector3d n = {v2.y() - v1.y(), v2.x() - v1.x(), 0.0};
    n.normalize();

    eps *= (v2 - v1).norm();

    for (const auto& v: nodes) {
        if (std::abs((v->pos() - c).dot(n)) > eps) {
            return false;
        }
    }
    return true;
}

void BlockStructured::set_boundary(const std::vector<BaseNode::Ptr>& nodes) {
    set_boundary(nodes, Boundary::UNDEFINED);
}

void BlockStructured::set_boundary(const std::vector<BaseNode::Ptr>& nodes, Boundary boundary) {
    if (m_stage != EDITABLE) {
        skip_warning("BlockStructured::set_boundary");
        return;
    }
    if (nodes.size() < 2) {
        throw std::invalid_argument("BlockStructured::set_boundary: at least two nodes required");
    }

    // Если образуют линию, то создаем её и устанавливаем
    if (collinear(nodes)) {
        Plane::Ptr plane = Plane::create(nodes.front(), nodes.back());
        if (boundary != Boundary::UNDEFINED) {
            plane->set_boundary(boundary);
        }
        set_boundary(nodes, plane);
        return;
    }

    for (int i = 1; i < nodes.size(); ++i) {
        set_boundary(nodes[i - 1], nodes[i], boundary);
    }
}

void BlockStructured::set_size(BaseNode::Ref v1, BaseNode::Ref v2, int size) {
    if (m_stage == EDITABLE) {
        link_blocks();
    }
    set_size(std::vector{v1, v2}, size);
}

void BlockStructured::set_size(const std::vector<BaseNode::Ptr>& nodes, int size) {
    if (m_stage == EDITABLE) {
        link_blocks();
    }
    if (nodes.size() < 2) {
        throw std::runtime_error("BlockStructured::set_size: can't set size, has no edges");
    }
    for (int i = 0; i < static_cast<int>(nodes.size() - 1); ++i) {
        BaseNode::Ref v1 = nodes[i];
        BaseNode::Ref v2 = nodes[i + 1];
        BaseEdge edge{v1, v2};
        if (!m_edges.contains(edge)) {
            std::cerr << "BlockStructured::set_size: can't set size, no such edge\n";
            return;
        }
    }
    m_wanted_sizes.push_back({.nodes=nodes, .size=size});
}

// ---------- Линковка ----------------------------------------- -------------------------------------------------------

void BlockStructured::remove_redundant() {
    const auto it = std::ranges::remove(m_blocks, nullptr).begin();
    m_blocks.erase(it, m_blocks.end());
    std::erase_if(m_edges, [](const auto& pair) {
        return pair.second.empty();
    });
    for (int i = 0; i < n_blocks(); ++i) {
        m_blocks[i]->set_index(i);
    }
}

// Проверка, что блоки связаны (в смысле граф связный)
inline bool connected_blocks(const std::vector<Block::Ptr>& blocks) {
    std::vector visited(blocks.size(), false);

    std::vector<int> to_visit;
    to_visit.reserve(blocks.size());
    to_visit.push_back(0);

    while (!to_visit.empty()) {
        int b1 = to_visit.back();
        to_visit.pop_back();
        visited[b1] = true;

        for (Side side: sides_2D) {
            auto neib = blocks[b1]->adjacent_block(side);
            if (neib) {
                int b2 = neib->index();
                if (!visited[b2]) {
                    to_visit.push_back(b2);
                }
            }
        }
    }

    bool all_visited = true;
    for (const auto v: visited) {
        if (!v) { all_visited = false; break;}
    }
    return all_visited;
}

void BlockStructured::link_blocks() {
    if (m_stage != EDITABLE) {
        throw std::runtime_error("BlockStructured::link: is only available at the editing stage");
    }

    // Удалить лишние элементы
    remove_redundant();

    // Собираем уникальные узлы
    std::map<BaseNode::Ptr, std::set<Block::Ptr>> nodes;
    for (const auto& block: m_blocks) {
        for (const auto& v: block->base_nodes()) {
            v->clear();
            nodes[v] = {};
        }
    }

    // Для каждого узла собираем смежные блоки
    for (auto& block: m_blocks) {
        for (const auto& v: block->base_nodes()) {
            nodes[v].insert(block);
        }
    }

    // Завершить формирование узлов (задать смежные блоки)
    for (auto& [node, adj]: nodes) {
        node->finalize(adj);
    }

    // Для каждой вершины проходим по смежным блокам
    for (const auto& node: nodes | std::views::keys) {
        if (node->n_adjacent_blocks() <= 1) {
            continue;
        }
        auto& adj = node->adjacent_blocks();
        for (int i = 0; i < adj.size(); ++i) {
            int j = (i + 1) % adj.size();
            Block::link(adj[i].lock(), adj[j].lock());
        }
    }

    // Проверка связности блоков
    if (!connected_blocks(m_blocks)) {
        throw std::runtime_error("BlockStructured::link_blocks() blocks is not connected");
    }

    // Оценить конформный модуль по геометрии
    for (const auto& block: m_blocks) {
        // При загрузке из кэша не меняем
        if (std::isnan(block->modulus())) {
            block->estimate_modulus();
        }
    }
    std::set<BaseNode::Ptr> inner_nodes_set;
    for (const auto& block: m_blocks) {
        for (const auto& node: block->base_nodes()) {
            // Пропускаем угловые и граничные базовые узлы
            if (node->fixed()) { continue; }

            if (node->n_adjacent_blocks() < 2) { continue; }

            if (auto [side1, side2] = block->incident_sides(node);
                block->boundary(side1) || block->boundary(side2)) {
                continue;
                }
            if (node->boundary()) {
                continue;
            }

            // Остались только внутренние базовые узлы
            inner_nodes_set.insert(node) = {};
        }
    }

    // Собрать внутренние узлы
    m_inner_nodes.clear();
    m_inner_nodes.reserve(inner_nodes_set.size());
    for (const auto& node: inner_nodes_set) {
        m_inner_nodes.emplace_back(node);
    }

    for (auto& [node, blocks]: m_inner_nodes) {
        // Предполагаем, что смежные вершины и блоки упорядочены верно,
        // то есть против часовой стрелки внутри области
        auto& adjacent_nodes  = node->adjacent_nodes();
        auto& adjacent_blocks = node->adjacent_blocks();

        if (adjacent_nodes.size() != adjacent_blocks.size()) {
            throw std::runtime_error("BlockStructured::consistent_modulus: n_blocks != n_nodes for inner node");
        }
        blocks.reserve(adjacent_blocks.size());
        for (int i = 0; i < adjacent_nodes.size(); ++i) {
            BaseNode::Ptr node2 = adjacent_nodes[i].lock();
            if (!node2) { throw std::runtime_error("BlockStructured::consistent_modulus: invalid node"); }
            const Block::Ptr block = adjacent_blocks[i].lock();
            blocks.emplace_back(block->index(), block->get_axis(node, node2) == Axis::X);
        }
    }

    m_stage = LINKED;
}

// ---------- Размеры блоков ----------------------------------------- -------------------------------------------------

inline void check_modulus(const std::vector<Block::Ptr>& blocks) {
    for (const auto& block: blocks) {
        if (bad_number(block->modulus())) {
            throw std::runtime_error("Need to define modulus before");
        }
    }
}

// Проставить размеры, в том числе у связанных блоков
template <typename Type>
void static_set_size(const std::vector<Block::Ptr>& blocks, Pairs<Type>& sizes, int b1, Axis axis, Type N) {
    z_assert(blocks.size() == sizes.size(), "static_set_size: different sizes");
    if (sizes[b1][axis] == N) { return; }
    sizes[b1][axis] = N;
    for (Side side: sides_by_axis(axis)) {
        auto neib = blocks[b1]->adjacent_block(side);
        if (neib) {
            int b2 = neib->index();
            auto twin = blocks[b1]->twin_face(side);
            if (sizes[b2][twin] != N) {
                static_set_size(blocks, sizes, b2, to_axis(twin), N);
            }
        }
    }
}

void BlockStructured::set_size(Pairs<int>& sizes, int b1, Axis axis, int N) const {
    static_set_size(m_blocks, sizes, b1, axis, N);
}

int BlockStructured::calc_cells(const Pairs<int>& sizes) const {
    int count = 0;
    for (int b1 = 0; b1 < n_blocks(); ++b1) {
        count += sizes[b1][Axis::X] * sizes[b1][Axis::Y];
    }
    return count;
}

int BlockStructured::calc_nodes(const Pairs<int>& sizes) const {
    int count = 0;
    for (int b1 = 0; b1 < n_blocks(); ++b1) {
        count += (sizes[b1][Axis::X] + 1) * (sizes[b1][Axis::Y] + 1);
    }
    return count;
}

int BlockStructured::calc_cells() const {
    if (m_blocks.empty()) {
        throw std::runtime_error("BlockStructured::calc_cells: has no blocks");
    }
    if (m_stage == EDITABLE) {
        // Должно быть выставлено автоматически в set_size();
        throw std::runtime_error("BlockStructured::calc_cells: can't calc cells on editing stage");
    }

    // Размеры по пожеланию пользователя
    auto sizes = wanted_block_sizes();
    return calc_cells(sizes);
}

int BlockStructured::calc_nodes() const {
    if (m_blocks.empty()) {
        throw std::runtime_error("BlockStructured::make: empty block structured");
    }
    if (m_stage == EDITABLE) {
        // Должно быть выставлено автоматически в set_size();
        throw std::runtime_error("BlockStructured::calc_nodes: can't calc nodes on editing stage");
    }

    // Размеры по пожеланию пользователя
    auto sizes = wanted_block_sizes();
    return calc_nodes(sizes);
}

void BlockStructured::set_rel_size(Pairs<double>& sizes, int b1, Axis axis, double N) const {
    static_set_size(m_blocks, sizes, b1, axis, N);
}

void BlockStructured::init_rel_sizes(Pairs<double>& sizes, int b1) const {
    if (bad_number(m_blocks[b1]->modulus())) {
        throw std::runtime_error("BlockStructured::init_rel_sizes: invalid modulus");
    }
    set_rel_size(sizes, b1, Axis::X, m_blocks[b1]->modulus());
    set_rel_size(sizes, b1, Axis::Y, 1.0);
}

bool BlockStructured::update_rel_sizes(Pairs<double>& sizes, int b1) const {
    if (bad_number(m_blocks[b1]->modulus())) {
        throw std::runtime_error("BlockStructured::update_rel_sizes: invalid modulus");
    }
    if (good_number(sizes[b1][Axis::X]) && good_number(sizes[b1][Axis::Y])) {
        return true;
    }
    if (bad_number(sizes[b1][Axis::X]) && good_number(sizes[b1][Axis::Y])) {
        set_rel_size(sizes, b1, Axis::X, sizes[b1][Axis::Y] * m_blocks[b1]->modulus());
    }
    if (bad_number(sizes[b1][Axis::Y]) && good_number(sizes[b1][Axis::X])) {
        set_rel_size(sizes, b1, Axis::Y, sizes[b1][Axis::X] / m_blocks[b1]->modulus());
    }
    return good_number(sizes[b1][Axis::X]) && good_number(sizes[b1][Axis::Y]);
}

void BlockStructured::propagate_sizes(Pairs<double>& rel_sizes) const {
    // Алгоритм не сработает для несвязных областей
    std::list<Block::Ptr> list;
    for (int i = 0; i < n_blocks(); ++i) {
        list.emplace_back(m_blocks[i]);
    }

    int count = 0;
    while (!list.empty() && count < 10 * n_blocks()) {
        const auto& block = list.front();
        list.pop_front();

        // Если размеры блока не определены полностью, то переместить в начало
        bool defined = update_rel_sizes(rel_sizes, block->index());
        if (!defined) {
            list.push_back(block);
        }
        ++count;
    }
    if (count > 10 * m_blocks.size()) {
        throw std::runtime_error("BlockStructured::propagate_sizes: infinite loop, probably not connected blocks");
    }
}

Pairs<int> BlockStructured::auto_block_sizes(int N) const {
    check_modulus(m_blocks);

    // Сбросить относительные размеры
    Pairs<double> rel_sizes(m_blocks.size(), {NAN, NAN});

    // Задаем у первого блока
    init_rel_sizes(rel_sizes, 0);

    // Распространить размеры на все блоки
    propagate_sizes(rel_sizes);

    // Находим минимальный относительный размер
    double min_size = std::numeric_limits<double>::max();
    for (int b1 = 0; b1 < n_blocks(); ++b1) {
        min_size = std::min(min_size, rel_sizes[b1][Axis::X]);
        min_size = std::min(min_size, rel_sizes[b1][Axis::Y]);
    }
    for (int b1 = 0; b1 < n_blocks(); ++b1) {
        rel_sizes[b1][Axis::X] *= (N / min_size);
        rel_sizes[b1][Axis::Y] *= (N / min_size);
    }

    // Проставить реальные размеры (N ячеек на min_size)
    Pairs<int> sizes(m_blocks.size());
    for (int b1 = 0; b1 < n_blocks(); ++b1) {
        sizes[b1][Axis::X] = static_cast<int>(std::round(rel_sizes[b1][Axis::X]));
        sizes[b1][Axis::Y] = static_cast<int>(std::round(rel_sizes[b1][Axis::Y]));
    }
    return sizes;
}

Pairs<int> BlockStructured::wanted_block_sizes() const {
    if (m_wanted_sizes.empty()) {
        throw std::runtime_error("BlockStructured::wanted_block_sizes: set at least one size");
    }

    // Все модули должны быть выставлены
    check_modulus(m_blocks);

    // Сбросить относительные размеры
    Pairs<double> rel_sizes(m_blocks.size(), {NAN, NAN});

    // Задаем размеры у части блоков
    for (const auto& [nodes, size]: m_wanted_sizes) {
        if (nodes.size() < 2) {
            throw std::runtime_error("BlockStructured::wanted_block_sizes: wrong node chain");
        }
        if (nodes.size() == 2) {
            BaseNode::Ptr v1 = nodes[0];
            BaseNode::Ptr v2 = nodes[1];
            BaseEdge edge{v1, v2};
            if (!m_edges.contains(edge)) {
                throw std::runtime_error("BlockStructured::wanted_block_sizes: wrong edge chain");
            }
            auto& pair = m_edges.find(edge)->second;
            if (pair.empty()) {
                throw std::runtime_error("BlockStructured::wanted_block_sizes: empty edge in chain");
            }
            set_rel_size(rel_sizes, pair.b1.lock()->index(), to_axis(pair.side1), size);
        }
        else {
            // Тут сложно, нужно правильно распределить число ячеек по всей цепочке
            // Число блоков, число ребер в цепочке
            int n = static_cast<int>(nodes.size()) - 1;

            // Собираем цепочку блоков
            std::vector<Block::Ptr> blocks(n);

            BaseEdge edge0{nodes[0], nodes[1]};
            if (!m_edges.contains(edge0)) {
                throw std::runtime_error("BlockStructured::wanted_block_sizes: wrong edge in chain");
            }
            blocks[0] = m_edges.find(edge0)->second.b1.lock();

            for (int i = 1; i < n; ++i) {
                Block::Ref prev = blocks[i - 1];
                BaseNode::Ptr v0 = nodes[i - 1];
                BaseNode::Ptr v1 = nodes[i];

                auto [w1, w2] = prev->adjacent_nodes(v1);
                if (w1 == v0) {
                    blocks[i] = prev->adjacent_block(prev->get_side(v1, w2));
                }
                else if (w2 == v0) {
                    blocks[i] = prev->adjacent_block(prev->get_side(v1, w1));
                }
                else {
                    throw std::runtime_error("BlockStructured::wanted_block_sizes: something wrong #1");
                }

                BaseNode::Ptr v2 = nodes[i + 1];
                BaseEdge edge{v1, v2};
                if (!m_edges.contains(edge)) {
                    throw std::runtime_error("BlockStructured::wanted_block_sizes: wrong edge in chain");
                }
                if (m_edges.find(edge)->second.b1.lock() != blocks[i]) {
                    throw std::runtime_error("BlockStructured::wanted_block_sizes: something wrong #2");
                }
            }

            std::vector<Axis> axes(n);
            for (int i = 0; i < n; ++i) {
                axes[i] = blocks[i]->get_axis(nodes[i], nodes[i + 1]);
            }

            // Небольшое число, локальные размеры
            Pairs<double> sizes(n);
            sizes[0][Axis::X] = blocks[0]->modulus();
            sizes[0][Axis::Y] = 1.0;
            for (int i = 1; i < n; ++i) {
                sizes[i] = sizes[i - 1];
                if (axes[i] != axes[i - 1]) {
                    std::swap(sizes[i][Axis::X], sizes[i][Axis::Y]);
                }
                if (axes[i] == Axis::X) {
                    sizes[i][Axis::X] = sizes[i][Axis::Y] * blocks[i]->modulus();
                }
                else {
                    sizes[i][Axis::Y] = sizes[i][Axis::X] / blocks[i]->modulus();
                }
            }

            double sum = 0.0;
            for (int i = 0; i < n; ++i) {
                sum += sizes[i][axes[i]];
            }

            for (int i = 0; i < n; ++i) {
                sizes[i][Axis::X] *= size / sum;
                sizes[i][Axis::Y] *= size / sum;
            }

            for (int i = 0; i < n; ++i) {
                int N = static_cast<int>(std::round(sizes[i][axes[i]]));
                set_rel_size(rel_sizes, blocks[i]->index(), axes[i], N);
            }
        }
    }

    // Распространить размеры на все блоки
    propagate_sizes(rel_sizes);

    // Целые размеры
    Pairs<int> sizes(m_blocks.size());
    for (int b1 = 0; b1 < n_blocks(); ++b1) {
        sizes[b1][Axis::X] = std::max(1, static_cast<int>(std::round(rel_sizes[b1][Axis::X])));
        sizes[b1][Axis::Y] = std::max(1, static_cast<int>(std::round(rel_sizes[b1][Axis::Y])));
    }

    return sizes;
}

void BlockStructured::sizes_info(const Pairs<int>& sizes) const {
    for (int b1 = 0; b1 < n_blocks(); ++b1) {
        std::cout << std::format("    Block {:3}.  K: {:.2f};  sizes: ({:3}, {:3})\n", b1,
            m_blocks[b1]->modulus(), sizes[b1][Axis::X], sizes[b1][Axis::Y]);
    }
}

void BlockStructured::conformal_info(const Pairs<double>& lambda) const {
    for (int b1 = 0; b1 < m_blocks.size(); ++b1) {
        std::cout << std::format("      Block {:3}.  K: {:.2f};  λ_1: {:.2f}; λ_2: {:.2f}\n",
            b1, m_blocks[b1]->modulus(), lambda[b1][Axis::X], lambda[b1][Axis::Y]);
    }
}

// ---------- Таблицы вершин, оптимизация, сглаживание, финальная сетка ------------------------------------------------

namespace {

// Проверить размеры смежных блоков перед созданием вершин
void check_consistency(const std::vector<Block::Ptr>& blocks, const Pairs<int>& sizes) {
    if (blocks.size() != sizes.size()) {
        throw std::runtime_error("BlockStructured::check_consistency: block sizes mismatch");
    }

    for (int b1 = 0; b1 < blocks.size(); ++b1) {
        const auto& block = blocks[b1];

        if (sizes[b1][Axis::X] < 1 || sizes[b1][Axis::Y] < 1) {
            throw std::runtime_error("BlockStructured::check_consistency: zero size of some block");
        }
        for (Side side: sides_2D) {
            if (block->boundary(side)) { continue; }

            auto neib = block->adjacent_block(side);
            if (neib) {
                int b2 = neib->index();
                Side twin = block->twin_face(side);

                int size1 = sizes[b1][side];
                int size2 = sizes[b2][twin];

                if (size1 != size2) {
                    throw std::runtime_error(std::format("Block::check_consistency: sizes mismatch (blocks {}, {})", b1, b2));
                }
            }
            else {
                throw std::runtime_error("Block::check_consistency: boundary or neighbor should be defined");
            }
        }
    }
}

// Сгенерировать таблицы узлов, если в блоках задано конформное отображение, то оно будет
// использовано для билинейной интерполяции. Затем связать вершины друг с другом.
Tables2D create_vertices(const std::vector<Block::Ptr>& blocks, const Pairs<int>& sizes, Pairs<double>& lambda) {
    check_consistency(blocks, sizes);

    // ---------- Создать таблицы вершин ----------------------------------------------------------
    Tables2D all_vertices(blocks.size());
    for (int b1 = 0; b1 < blocks.size(); ++b1) {
        all_vertices[b1] = blocks[b1]->create_vertices(sizes[b1]);
    }

    // ---------- Склеить вершины на границах блоков-----------------------------------------------
    for (int b1 = 0; b1 < blocks.size(); ++b1) {
        const auto& block = blocks[b1];

        // Связать вершины с соседом
        for (Side side: sides_2D) {
            auto neib = block->adjacent_block(side);

            if (!neib) { continue; }

            int b2 = neib->index();

            // Сосед есть, но без вершин
            if (all_vertices[b2].empty()) {
                continue;
            }

            Side twin = block->twin_face(side);

            int N = sizes[b1][side];
            for (int idx = 0; idx <= N; ++idx) {
                all_vertices[b1].boundary(side, idx) = all_vertices[b2].boundary(twin, N - idx);
            }
        }
    }

    // ---------- Склеить угловые вершины блоков---------------------------------------------------
    std::set<BaseNode::Ptr> unique_nodes;
    for (const auto & block : blocks) {
        unique_nodes.insert(block->base_node(0));
        unique_nodes.insert(block->base_node(1));
        unique_nodes.insert(block->base_node(2));
        unique_nodes.insert(block->base_node(3));
    }
    for (const auto& node: unique_nodes) {
        auto& adj_blocks = node->adjacent_blocks();
        if (adj_blocks.size() < 2) { continue; }

        Block::Ptr base_block = adj_blocks.front().lock();
        int base_v_idx = base_block->base_node_index(node);
        BsVertex::Ptr main_corner = all_vertices[base_block->index()].corner(base_v_idx);
        for (int i = 1; i < adj_blocks.size(); ++i) {
            Block::Ptr neib_block = adj_blocks[i].lock();
            int neib_v_idx = neib_block->base_node_index(node);
            all_vertices[neib_block->index()].corner(neib_v_idx) = main_corner;
        }
    }

    // ---------- Связать вершины друг с другом ----------------------------------------------------
    for (int b1 = 0; b1 < blocks.size(); ++b1) {
        const auto& self = blocks[b1];
        auto& vertices = all_vertices[b1];

        // Очистить все списки
        for (int i = 0; i <= sizes[b1][Axis::X]; ++i) {
            for (int j = 0; j <= sizes[b1][Axis::Y]; ++j) {
                vertices(i, j)->clear();
            }
        }

        // Внутренние вершины блоков
        for (int i = 1; i < sizes[b1][Axis::X]; ++i) {
            for (int j = 1; j < sizes[b1][Axis::Y]; ++j) {
                // Справа, сверху, слева, снизу
                std::vector edges {
                    BsEdge::Inside(vertices(i + 1, j), &lambda[b1][Axis::X]),
                    BsEdge::Inside(vertices(i, j + 1), &lambda[b1][Axis::Y]),
                    BsEdge::Inside(vertices(i - 1, j), &lambda[b1][Axis::X]),
                    BsEdge::Inside(vertices(i, j - 1), &lambda[b1][Axis::Y])
                };
                vertices(i, j)->set_edges(edges);
            }
        }

        // Вершины на границах блока (без угловых)
        for (Side side: sides_2D) {
            if (!self->boundary(side) && !self->adjacent_block(side)) {
                throw std::runtime_error("create_vertices: face is not boundary and has no neighbor #1");
            }

            Axis axis = to_axis(side);
            Axis perp_axis = orthogonal(axis);

            if (self->boundary(side)) {
                for (int idx = 1; idx < sizes[b1][side]; ++idx) {
                    // Добавить границу
                    vertices.boundary(side, idx)->add_boundary(self->boundary(side).get());
                    // Добавить связи (три штуки)
                    std::vector edges {
                        BsEdge::Border(vertices.boundary(side, idx + 1),  &lambda[b1][axis]),
                        BsEdge::Inside(vertices.near_boundary(side, idx), &lambda[b1][perp_axis]),
                        BsEdge::Border(vertices.boundary(side, idx - 1),  &lambda[b1][axis])
                    };
                    vertices.boundary(side, idx)->set_edges(edges);
                }
            } else {
                auto neib = self->adjacent_block(side);
                if (!neib) {
                    throw std::runtime_error("create_vertices: face is not boundary and has no neighbor #2");
                }
                int b2 = neib->index();

                int N = sizes[b1][side];
                Side twin = self->twin_face(side);

                Axis neib_axis = to_axis(twin);
                Axis perp_neib_axis = orthogonal(neib_axis);
                for (int idx = 1; idx < N; ++idx) {
                    std::vector edges {
                        BsEdge::Inside(vertices.boundary(side, idx + 1),  &lambda[b1][axis], &lambda[b2][neib_axis]),
                        BsEdge::Inside(vertices.near_boundary(side, idx), &lambda[b1][perp_axis]),
                        BsEdge::Inside(vertices.boundary(side, idx - 1),  &lambda[b1][axis], &lambda[b2][neib_axis]),
                        BsEdge::Inside(all_vertices[b2].near_boundary(twin, N - idx), &lambda[b2][perp_neib_axis])
                    };
                    vertices.boundary(side, idx)->set_edges(edges);
                }
            }
        }

        // Угловые вершины блока (здесь могут быть сингулярности)
        for (int v_idx = 0; v_idx < 4; ++v_idx) {
            auto node = self->base_node(v_idx);
            vertices.corner(v_idx)->clear();

            /// Фиксированная точка
            if (node->fixed()) { continue; }

            // Угловая вершина области, должна быть неподвижной
            if (node->n_adjacent_blocks() < 2) { continue; }

            // Угловая вершина области по другому критерию
            if (auto [side1, side2] = node_incident_sides(v_idx);
                self->boundary(side1) && self->boundary(side2)) {
                continue;
                }

            // Вершина на границе
            if (node->boundary()) {
                // Предполагаем, что смежные вершины и блоки упорядочены верно,
                // то есть против часовой стрелки внутри области
                auto& adjacent_nodes = node->adjacent_nodes();
                auto& adjacent_blocks = node->adjacent_blocks();
                if (adjacent_blocks.size() + 1 != adjacent_nodes.size()) {
                    throw std::runtime_error("create_vertices: n_blocks + 1 != n_nodes for boundary node");
                }

                auto corner = vertices.corner(v_idx);
                std::vector<BsEdge> edges;

                // Первая граничная вершина
                {
                    BaseNode::Ptr node2 = adjacent_nodes.front().lock();
                    Block::Ptr block = adjacent_blocks.front().lock();
                    if (!node2 || !block) { throw std::runtime_error("create_vertices: invalid block"); }
                    auto side = block->get_side(node, node2);
                    if (!block->boundary(side)) {
                        throw std::runtime_error("create_vertices: invalid block (not boundary)");
                    }
                    corner->add_boundary(block->boundary(side).get());
                    edges.push_back(BsEdge::Border(
                        all_vertices[block->index()].boundary(side, 1),
                        &lambda[block->index()][side]
                    ));
                }
                // Внутренние вершины
                for (int i = 1; i < adjacent_blocks.size(); ++i) {
                    BaseNode::Ptr node2 = adjacent_nodes[i].lock();
                    if (!node2) { throw std::runtime_error("create_vertices: invalid node"); }
                    Block::Ptr block = adjacent_blocks[i].lock();
                    auto prev_block = adjacent_blocks[i - 1].lock();
                    if (!block || !prev_block) { throw std::runtime_error("create_vertices: invalid block"); }
                    int b3 = block->index();
                    Side side = block->get_side(node, node2);
                    Side prev_side = prev_block->get_side(node, node2);
                    edges.push_back(BsEdge::Inside(
                        all_vertices[b3].boundary(side, 1),
                        &lambda[b3][side],
                        &lambda[prev_block->index()][prev_side]
                    ));
                    if (block->boundary(side)) {
                        throw std::runtime_error("create_vertices: invalid block (why boundary)");
                    }
                }
                // Последняя граничная вершина
                {
                    BaseNode::Ptr node2 = adjacent_nodes.back().lock();
                    Block::Ptr block = adjacent_blocks.back().lock();
                    if (!node2 || !block) { throw std::runtime_error("create_vertices: invalid block"); }
                    Side side = block->get_side(node, node2);
                    edges.push_back(BsEdge::Border(
                        all_vertices[block->index()].boundary(side, sizes[block->index()][side] - 1),
                        &lambda[block->index()][side]
                    ));
                    if (!block->boundary(side)) {
                        throw std::runtime_error("create_vertices: invalid block (not boundary)");
                    }
                    corner->add_boundary(block->boundary(side).get());
                }
                corner->set_edges(edges);
            }
            else {
                // Предполагаем, что смежные вершины и блоки упорядочены верно,
                // то есть против часовой стрелки внутри области
                auto& adjacent_nodes = node->adjacent_nodes();
                auto& adjacent_blocks = node->adjacent_blocks();

                if (adjacent_nodes.size() != adjacent_blocks.size()) {
                    throw std::runtime_error("create_vertices: n_blocks != n_nodes for inner node");
                }

                auto corner = vertices.corner(v_idx);
                std::vector<BsEdge> edges;

                for (int i = 0; i < adjacent_blocks.size(); ++i) {
                    BaseNode::Ptr node2 = adjacent_nodes[i].lock();
                    if (!node2) { throw std::runtime_error("create_vertices: invalid node"); }
                    Block::Ptr block = adjacent_blocks[i].lock();
                    int j = int(i + adjacent_blocks.size() - 1) % adjacent_blocks.size();
                    auto prev_block = adjacent_blocks[j].lock();
                    if (!block || !prev_block) { throw std::runtime_error("create_vertices: invalid block"); }
                    Side side = block->get_side(node, node2);
                    Side prev_side = prev_block->get_side(node, node2);
                    edges.push_back(BsEdge::Inside(
                        all_vertices[block->index()].boundary(side, 1),
                        &lambda[block->index()][side],
                        &lambda[prev_block->index()][prev_side]
                    ));
                    if (block->boundary(side)) {
                        throw std::runtime_error("create_vertices: invalid block (why boundary)");
                    }
                }
                corner->set_edges(edges);
            }
        }
    }

    return all_vertices;
}

// Сглаживание вершин в блоке
void smooth_one(const Table2D& vertices) {
    for (int i = 0; i < vertices.size(Axis::X); ++i) {
        for (int j = 0; j < vertices.size(Axis::Y); ++j) {
            const auto& vc = vertices(i, j);
            auto &adjacent = vc->adjacent();

            // Фиксированная вершина
            if (adjacent.empty()) {
                vc->next =vc->pos;
                continue;
            }

            // Угловая вершина
            if (!vc->inner() && vc->corner()) {
                vc->next = vc->pos;
                continue;
            }

            // Вершина на границе
            if (!vc->inner()) {
                Curve* curve = vc->boundary();
                Vector3d next = Vector3d::Zero();
                double sum = 0.0;
                for (auto edge: vc->adjacent()) {
                    if (edge.boundary()) {
                        next += edge.lambda() * edge.pos();
                        sum += edge.lambda();
                    }
                    else {
                        next += 2.0 * edge.lambda() * curve->projection(edge.pos());
                        sum += 2.0 * edge.lambda();
                    }
                }
                next /= sum;
                vc->next = curve->projection(next);
                continue;
            }

            // Внутренняя вершина
            Vector3d next = Vector3d::Zero();
            double sum = 0.0;
            for (auto edge: vc->adjacent()) {
                next += edge.lambda() * edge.pos();
                sum += edge.lambda();
            }
            next /= sum;
            vc->next = next;
        }
    }
}

// Сглаживание вершин во всех блоках блоке.
// Возвращает максимальный относительный сдвиг вершин.
double smooth_vertices(const Tables2D& all_vertices) {
    for (const auto& vertices: all_vertices) {
        smooth_one(vertices);
    }

    double err = 0.0;
    for (auto& vertices: all_vertices) {
        for (int i = 0; i < vertices.size(Axis::X); ++i) {
            for (int j = 0; j < vertices.size(Axis::Y); ++j) {
                BsVertex::Ref vc = vertices(i, j);
                double L = 0.0;
                for (const auto edge: vc->adjacent()) {
                    L = std::max(L, (vc->pos - edge.pos()).norm());
                }
                err = std::max(err, (vc->pos - vc->next).norm() / L);
                vc->pos = vc->next;
            }
        }
    }
    return err;
}

} // anonymous namespace

inline AxisPair<double> calc_lambda(AxisPair<int> sizes, double K) {
    AxisPair<double> lambda{NAN, NAN};
    if (sizes[Axis::X] >= 1 && sizes[Axis::Y] >= 1) {
        lambda[Axis::Y] = K * sizes[Axis::Y] / sizes[Axis::X];
        lambda[Axis::X] = 1.0 / lambda[Axis::Y];
    }
    return lambda;
}

// Коэффициенты сглаживания (на основе размеров и конформных модулей)
inline Pairs<double> get_lambda(const std::vector<Block::Ptr>& blocks, const Pairs<int>& sizes) {
    if (blocks.size() != sizes.size()) {
        throw std::runtime_error("get_lambda: block sizes mismatch");
    }
    for (int b1 = 0; b1 < blocks.size(); ++b1) {
    }
    Pairs<double> lambda(blocks.size());
    for (int b1 = 0; b1 < blocks.size(); ++b1) {
        if (bad_number(blocks[b1]->modulus()) || sizes[b1][Axis::X] < 1 || sizes[b1][Axis::Y] < 1) {
            throw std::runtime_error("get_lambda: bad block size");
        }
        lambda[b1] = calc_lambda(sizes[b1], blocks[b1]->modulus());
    }
    return lambda;
}

const std::string cachefile = ".bscache";

template <class T>
inline void hash_combine(std::size_t& seed, const T& v) {
    seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

size_t BlockStructured::get_hash(const optimize_options& opts) const {
    // Считаем, что у нас топология не меняется, все кривые и границы заданы так же.
    // В контрольную сумму попадают только координаты базовых вершин и статус fix/не fix.
    size_t hash = m_blocks.size();

    // Хэшируем опции, они влияют на итоговую сетку
    hash_combine(hash, opts.steps);
    hash_combine(hash, opts.N);
    hash_combine(hash, opts.eps);

    for (const auto& block: m_blocks) {
        for (const auto& node: block->base_nodes()) {
            hash_combine(hash, node->x());
            hash_combine(hash, node->y());
            hash_combine(hash, node->fixed());
        }
    }
    return hash;
}

void BlockStructured::save_cache(const optimize_options& opts) const {
    size_t hash = get_hash(opts);
    std::ofstream file(cachefile, std::ios::out | std::ios::binary | std::ios::trunc);
    file.write((char*)&hash, sizeof(hash));

    for (const auto& block: m_blocks) {
        BaseNode::Ptr v1 = block->base_node(0);
        BaseNode::Ptr v2 = block->base_node(1);
        BaseNode::Ptr v3 = block->base_node(2);
        BaseNode::Ptr v4 = block->base_node(3);

        // Положения базисных узлов
        file.write((char*)v1->pos().data(), 3 * sizeof(double));
        file.write((char*)v2->pos().data(), 3 * sizeof(double));
        file.write((char*)v3->pos().data(), 3 * sizeof(double));
        file.write((char*)v4->pos().data(), 3 * sizeof(double));

        // Конформный модуль
        double modulus = block->modulus();
        file.write((char*)&modulus, sizeof(modulus));

        // Буфер вершин mapping
        std::array<int, 2> sizes = block->mapping().sizes();
        file.write((char*)sizes.data(), sizeof(sizes));
        auto buffer = block->mapping().data();
        file.write((char*)buffer, 3 * sizes[0] * sizes[1] * sizeof(double));
    }
    file.close();
}

void BlockStructured::load_cached() {
    std::ifstream file(cachefile, std::ios::in | std::ios::binary);
    size_t hash;
    file.read((char*)&hash, sizeof(hash));

    for (const auto& block: m_blocks) {
        // Положения базисных узлов
        Vector3d v1, v2, v3, v4;
        file.read((char*)v1.data(), 3 * sizeof(double));
        file.read((char*)v2.data(), 3 * sizeof(double));
        file.read((char*)v3.data(), 3 * sizeof(double));
        file.read((char*)v4.data(), 3 * sizeof(double));

        block->base_node(0)->set_pos(v1);
        block->base_node(1)->set_pos(v2);
        block->base_node(2)->set_pos(v3);
        block->base_node(3)->set_pos(v4);

        // Конформный модуль
        double modulus;
        file.read((char*)&modulus, sizeof(modulus));
        block->set_modulus(modulus);

        // Буфер вершин mapping
        std::array<int, 2> sizes;
        file.read((char*)sizes.data(), sizeof(sizes));

        std::vector<Vector3d> buffer(sizes[0]*sizes[1]);
        file.read((char*)buffer.data(), 3 * buffer.size() * sizeof(double));
        Array2D<Vector3d> vertices(sizes, std::move(buffer));
        block->set_mapping(std::move(vertices));
    }
    file.close();

    if (m_stage == EDITABLE) {
        if (m_verbosity > 1) {
            std::cout << "  Link blocks\n";
        }
        link_blocks();
    }
    m_stage = OPTIMIZED;
}

void BlockStructured::optimize(const optimize_options& opts) {
    namespace fs = std::filesystem;

    int n_steps = std::max(1, std::min(opts.steps, 7));
    int N = std::max(2, opts.N);

    // На последней итерации N не должно превышать 1024
    N = std::min(N, static_cast<int>(std::floor(std::pow(2, 11 - n_steps))));

    double eps = std::max(1.0e-6, std::min(opts.eps, 0.01));

    // Загрузить из кэша, если там что-то есть
    if (opts.cache) {
        // Проверить наличие файла кэша
        if (fs::exists(cachefile) && fs::is_regular_file(cachefile)) {
            std::ifstream file(cachefile, std::ios::in | std::ios::binary);
            size_t file_hash;
            file.read((char*)&file_hash, sizeof(file_hash));
            file.close();
            size_t hash = get_hash(opts);
            if (file_hash == hash) {
                if (m_verbosity > 0) {
                    std::cout << "BlockStructured: load cached \"" << cachefile << "\"\n";
                }
                load_cached();
                return;
            }
        }
    }

    if (m_verbosity > 0) {
        std::cout << "BlockStructured optimization (" << n_steps << " steps)\n";
        std::cout << "  verbosity: " << m_verbosity << "\n";
        std::cout << "  epsilon:   " << std::format("{:.2e};\n", eps);
        std::cout << "  start from " << N << " cells.\n";
    }
    for (int step = 0; step < n_steps; ++step) {
        if (m_verbosity > 0) {
            std::cout << std::format("Optimization step {}.  N: {}\n", step + 1, N);
        }
        optimize(N, eps);
        N *= 2;
        eps /= 2.0;
    }

    // Сохранить кэш после оптимизации
    if (opts.cache) {
        if (m_verbosity > 0) {
            std::cout << "BlockStructured: save cache \"" << cachefile << "\"\n";
        }
        save_cache(opts);
    }
}

void BlockStructured::optimize(int N, double eps) {
    if (m_blocks.empty()) {
        throw std::runtime_error("BlockStructured::optimize: add blocks before optimization");
    }
    if (m_stage == EDITABLE) {
        if (m_verbosity > 1) {
            std::cout << "  Link blocks\n";
        }
        link_blocks();
    }

    if (m_verbosity > 1) {
        if (m_stage == LINKED) {
            std::cout << "  Init optimization\n";
        } else {
            std::cout << "  Repeat optimization\n";
        }
    }
    optimization_step(N, eps);
}

void BlockStructured::optimization_step(int N, double eps) {
    // Проставить размеры блоков
    const auto sizes = auto_block_sizes(N);
    auto lambda = get_lambda(m_blocks, sizes);

    if (m_verbosity > 1) {
        std::cout << "    Cells count: " << calc_cells(sizes) << "\n";
        if (m_verbosity > 2) {
            sizes_info(sizes);
        }
    }

    // Создание вершин
    const auto vertices = create_vertices(m_blocks, sizes, lambda);
    if (m_verbosity > 4 && m_blocks[0]->mapping().empty()) {
        // Показать при первом создании
        plot(vertices);
    }

    int counter = 0;
    double error = 1.0;

    while (error > eps && counter < 5000) {
        if (m_verbosity > 2 && counter % 50 == 0) {
            std::cout << std::format("    step {:6}, eps: {:.2e}\n", counter, error);
            if (m_verbosity > 3 && counter % 50 == 0) {
                conformal_info(lambda);
            }
        }
        error = smooth_vertices(vertices);
        for (int b1 = 0; b1 < n_blocks(); ++b1) {
            m_blocks[b1]->update_modulus(vertices[b1]);
        }
        correct_modulus();
        for (int b1 = 0; b1 < n_blocks(); ++b1) {
            lambda[b1] = calc_lambda(sizes[b1], m_blocks[b1]->modulus());
        }
        ++counter;
    }
    if (m_verbosity > 2) {
        std::cout << std::format("    step {:6}, eps: {:.2e}\n", counter, error);
        if (m_verbosity > 3) {
            conformal_info(lambda);
        }
    }
    if (m_verbosity > 4) {
        plot(vertices);
    }
    for (int k = 0; k < n_blocks(); ++k) {
        m_blocks[k]->set_mapping(vertices[k]);
    }
    m_stage = OPTIMIZED;
}

void BlockStructured::correct_modulus() {
    using Eigen::VectorXd;
    using Eigen::MatrixXd;

    // Решается недоопределенная система уравнений
    // Число строк - число внутренних узлов, столбцов - число блоков.
    int n_rows = static_cast<int>(m_inner_nodes.size());

    MatrixXd A = MatrixXd::Zero(n_rows, n_blocks());
    VectorXd F = VectorXd::Zero(n_rows);

    for (int i = 0; i < n_rows; ++i) {
        for (auto [b, dir]: m_inner_nodes[i].blocks) {
            double sign = dir ? 1.0 : -1.0;
            A(i, b) = sign;
            F[i] -= sign * std::log(m_blocks[b]->modulus());
        }
    }

    // Псевдообратная матрица
    MatrixXd C = A * A.transpose();
    VectorXd y = C.ldlt().solve(F);
    VectorXd x = A.transpose() * y;

    VectorXd xi(n_blocks());
    for (int b1 = 0; b1 < m_blocks.size(); ++b1) {
        xi[b1] = std::exp(x[b1]);
    }

    for (int b1 = 0; b1 < n_blocks(); ++b1) {
        m_blocks[b1]->set_modulus(xi[b1] * m_blocks[b1]->modulus());
    }

    /*
    // Debug info
    for (auto& [node, info]: m_inner_nodes) {
        std::cout << "  Node (" << std::boolalpha;
        for (auto [b1, dir]: info) {
            std::cout << b1 << ",";
        }
        std::cout << ")\t";
        double prod = 1.0;
        for (auto [b1, dir]: info) {
            double sign = dir ? 1.0 : -1.0;
            prod *= std::pow(m_blocks[b1]->modulus(), sign);
        }
        std::cout << "prod: " << prod << "\n";
    }
    */
}

// Пронумеровать вершины в блоках, возвращает число уникальных вершин
inline int enumerate(const Tables2D& vertices) {
    for (const auto& verts: vertices) {
        for (int i = 0; i < verts.size(Axis::X); ++i) {
            for (int j = 0; j < verts.size(Axis::Y); ++j) {
                verts(i, j)->index = -1;
            }
        }
    }

    int n_nodes = 0;
    for (const auto& verts: vertices) {
        for (int i = 0; i < verts.size(Axis::X); ++i) {
            for (int j = 0; j < verts.size(Axis::Y); ++j) {
                if (verts(i, j)->index < 0) {
                    verts(i, j)->index = n_nodes;
                    ++n_nodes;
                }
            }
        }
    }
    return n_nodes;
}

struct Smoother {
    int n_nodes;  // Число узлов

    std::vector<Vector3d> pos;  // Текущие положения вершин
    std::vector<Vector3d> next; // Следующие положения вершин

    Csr<int, int>    neibs;   // Индексы соседних вершин
    Csr<double, int> lambda;  // Коэффициенты сглаживания

    std::vector<Curve*> curves;    // Указатели на границы
    Csr<int, int> boundary_nodes;  // Индексы узлов, принадлежащих границам

    Smoother(const Tables2D& vertices, int nodes_count) : n_nodes(nodes_count) {
        pos.resize(n_nodes);
        std::vector<int> node_degrees(n_nodes);
        std::unordered_map<Curve*, std::vector<int>> boundary_indices;
        for (const auto& verts: vertices) {
            for (int i = 0; i < verts.size(Axis::X); ++i) {
                for (int j = 0; j < verts.size(Axis::Y); ++j) {
                    const auto& vc = verts(i, j);
                    int node_idx = vc->index;
                    pos[node_idx] = vc->pos;

                    if (vc->corner()) {
                        node_degrees[node_idx] = 0;
                    }
                    else {
                        node_degrees[node_idx] = vc->degree();
                        if (vc->is_boundary()) {
                            boundary_indices[vc->boundary()].push_back(node_idx);
                        }
                    }
                }
            }
        }

        // Скопируем, угловые вообще зафиксированы
        next = pos;

        // Расширим массивы под соседей
        neibs .assign_row_sizes(node_degrees);
        lambda.assign_row_sizes(node_degrees);

        auto restricted_lambda = [](double lambda) -> double {
            constexpr double threshold = 1.1;
            if (lambda < threshold || lambda > 1.0 / threshold) {
                return 1.0;
            }
            return lambda;
        };

        for (const auto& verts: vertices) {
            for (int i = 0; i < verts.size(Axis::X); ++i) {
                for (int j = 0; j < verts.size(Axis::Y); ++j) {
                    const auto& vc = verts(i, j);

                    if (vc->corner()) { continue; }

                    int node_idx = vc->index;

                    auto& edges = vc->adjacent();
                    double sum = 0.0;
                    for (int k = 0; k < edges.size(); ++k) {
                        neibs [node_idx][k] = edges[k].neib_idx();
                        sum += restricted_lambda(edges[k].lambda());
                    }
                    for (int k = 0; k < edges.size(); ++k) {
                        lambda[node_idx][k] = restricted_lambda(edges[k].lambda()) / sum;
                    }
                }
            }
        }

        // Число вершин на каждую кривую
        std::vector<int> n_nodes_for_curve;
        for (const auto& [curve, ids]: boundary_indices) {
            curves.emplace_back(curve);
            n_nodes_for_curve.emplace_back(ids.size());
        }

        boundary_nodes.assign_row_sizes(n_nodes_for_curve);

        for (int c = 0; c < curves.size(); ++c) {
            const auto& ids = boundary_indices[curves[c]];
            for (int i = 0; i < n_nodes_for_curve[c]; ++i) {
                boundary_nodes[c][i] = ids[i];
            }
        }
    }

    void smooth() {
        utils::threads::parallel_for(0, n_nodes,
            [this](int i) {
                if (neibs[i].empty()) { return; }
                next[i] = Vector3d::Zero();
                for (int j = 0; j < neibs[i].size(); ++j) {
                    next[i] += lambda[i][j] * pos[neibs[i][j]];
                }
            });

        for (int c = 0; c < curves.size(); ++c) {
            Curve* curve = curves[c];
            for (int i: boundary_nodes[c]) {
                next[i] = curve->projection(next[i]);
            }
        }
        std::swap(pos, next);
    }

    // Очищает всё кроме pos
    void clean() {
        next.clear();
        neibs.clear();
        lambda.clear();
        curves.clear();
        boundary_nodes.clear();
        next.shrink_to_fit();
        neibs.shrink_to_fit();
        lambda.shrink_to_fit();
        curves.shrink_to_fit();
        boundary_nodes.shrink_to_fit();
    }
};

Grid BlockStructured::make() const {
    if (m_blocks.empty()) {
        throw std::runtime_error("BlockStructured::make: empty block structured");
    }
    if (m_stage == EDITABLE) {
        // Должно быть выставлено автоматически в set_size();
        throw std::runtime_error("BlockStructured::make: can't make block structured before link_blocks");
    }
    if (m_verbosity > 0) {
        std::cout << "BlockStructured make grid\n";
    }

    // Размеры по пожеланию пользователя
    auto sizes = wanted_block_sizes();

    // Нелинейная адаптивная сетка (в 4 раза больше ячеек)
    if (m_adaptive && !m_linear) {
        for (auto& size: sizes) {
            size[Axis::X] *= 2;
            size[Axis::Y] *= 2;
        }
    }

    // Рассчитать коэффициенты
    auto lambda = get_lambda(m_blocks, sizes);
    const auto vertices = create_vertices(m_blocks, sizes, lambda);

    int n_cells = calc_cells(sizes);
    int n_nodes = enumerate(vertices);

    if (m_verbosity > 1) {
        std::cout << "  Cells count: " << n_cells << "\n";
        std::cout << "  Nodes count: " << n_nodes << "\n";
        if (m_verbosity > 2) {
            sizes_info(sizes);
        }
    }

    std::vector<Node::Ptr> nodes;
    {
        if (m_verbosity > 0) {
            std::cout << "  Smooth final grid (" << m_iters_count << " iterations)\n";
        }
        Smoother s(vertices, n_nodes);
        for (int counter = 0; counter < m_iters_count; ++counter) {
            s.smooth();
            if (m_verbosity > 2 && counter % 500 == 0) {
                std::cout << std::format("    step {:6}\n", counter);
            }
        }
        s.clean();

        nodes.resize(n_nodes);
        for (int i = 0; i < n_nodes; ++i) {
            nodes[i] = Node::create(s.pos[i]);
        }
    }

    Grid grid;
    grid.reserve_nodes(n_nodes);
    grid.reserve_cells(n_cells);

    for (int i = 0; i < n_nodes; ++i) {
        grid.add_node(nodes[i]);
    }

    std::vector block_bc = {
        Boundary::INNER,
        Boundary::INNER,
        Boundary::INNER,
        Boundary::INNER,
    };
    std::vector bc = block_bc;

    for (int b1 = 0; b1 < n_blocks(); ++b1) {
        auto block = m_blocks[b1];
        for (Side side: sides_2D) {
            if (block->boundary(side)) {
                block_bc[static_cast<int>(side)] = block->boundary(side)->boundary();
            }
            else {
                block_bc[static_cast<int>(side)] = Boundary::INNER;
            }
        }

        const auto& verts = vertices[b1];
        // Не адаптивная или адаптивная линейная
        if (!m_adaptive || (m_adaptive && m_linear)) {
            int Nx = verts.size(Axis::X) - 1;
            int Ny = verts.size(Axis::Y) - 1;

            // Обход против часовой стрелки, с нижней левой вершины и нижней грани
            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    bc[0] = bc[1] = bc[2] = bc[3] = Boundary::INNER;
                    if (i == 0)      bc[Side2D::L] = block_bc[static_cast<int>(Side::L)];
                    if (i == Nx - 1) bc[Side2D::R] = block_bc[static_cast<int>(Side::R)];
                    if (j == 0)      bc[Side2D::B] = block_bc[static_cast<int>(Side::B)];
                    if (j == Ny - 1) bc[Side2D::T] = block_bc[static_cast<int>(Side::T)];

                    grid.add_cell(
                        CellType::QUAD, {
                            nodes[verts(i, j)->index],
                            nodes[verts(i + 1, j)->index],
                            nodes[verts(i + 1, j + 1)->index],
                            nodes[verts(i, j + 1)->index]
                        }, bc);
                }
            }
        }
        else if (m_adaptive && !m_linear) {
            // Адаптивная и нелинейная
            // AMR-ячейка, Z-порядок вершин
            int Nx = (verts.size(Axis::X) - 1) / 2;
            int Ny = (verts.size(Axis::Y) - 1) / 2;

            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    bc[0] = bc[1] = bc[2] = bc[3] = Boundary::INNER;
                    if (i == 0)      bc[Side2D::L] = block_bc[static_cast<int>(Side::L)];
                    if (i == Nx - 1) bc[Side2D::R] = block_bc[static_cast<int>(Side::R)];
                    if (j == 0)      bc[Side2D::B] = block_bc[static_cast<int>(Side::B)];
                    if (j == Ny - 1) bc[Side2D::T] = block_bc[static_cast<int>(Side::T)];

                    grid.add_cell(
                        CellType::AMR2D, {
                            nodes[verts(2*i, 2*j+0)->index], nodes[verts(2*i+1, 2*j+0)->index], nodes[verts(2*i+2, 2*j+0)->index],
                            nodes[verts(2*i, 2*j+1)->index], nodes[verts(2*i+1, 2*j+1)->index], nodes[verts(2*i+2, 2*j+1)->index],
                            nodes[verts(2*i, 2*j+2)->index], nodes[verts(2*i+1, 2*j+2)->index], nodes[verts(2*i+2, 2*j+2)->index]
                        }, bc);
                }
            }
        }
        else {
            throw std::runtime_error("BlockStructured::make: bad configuration");
        }
    }

    // Сделать адаптивную из простой сетки
    if (m_adaptive && m_linear) {
        grid.make_amr();
    }

    if (m_verbosity > 0) {
        std::cout << "BlockStructured finish\n";
    }
    return grid;
}

} // namespace zephyr::geom::generator