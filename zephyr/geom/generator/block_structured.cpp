#include <array>
#include <algorithm>
#include <format>
#include <map>
#include <iomanip>
#include <list>
#include <ranges>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>
#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/block_structured.h>
#include <zephyr/utils/pyplot.h>

namespace zephyr::geom::generator {

inline bool good_number(double num) {
    return !std::isinf(num) && !std::isnan(num) && num > 0.0;
}

inline bool bad_number(double num) {
    return std::isinf(num) || std::isnan(num) || num <= 0.0;
}

constexpr int node_idx(Side side, int idx) {
    switch (side) {
        case Side::L: return std::array{0, 2}[idx];
        case Side::R: return std::array{1, 3}[idx];
        case Side::B: return std::array{0, 1}[idx];
        case Side::T: return std::array{2, 3}[idx];
        default: throw std::runtime_error("Invalid side");
    }
}

// Перпендикулярная ось
constexpr Axis orthogonal(Axis axis) {
    return axis == Axis::X ? Axis::Y : Axis::X;
}

// Противоположная сторона
constexpr Side opposite(Side side) {
    switch (side) {
        case Side::L: return Side::R;
        case Side::R: return Side::L;
        case Side::B: return Side::T;
        case Side::T: return Side::B;
        default: throw std::runtime_error("Invalid side");
    }
}

// Стороны, которые лежат вдоль выбранной оси
constexpr std::array<Side, 2> sides_by_axis(Axis axis) {
    return axis == Axis::X? std::array{Side::B, Side::T} : std::array{Side::L, Side::R};
}

// Пара сторон, связанных с вершиной (обход против часовой внутри блока)
constexpr std::tuple<Side, Side> node_adjacent_sides(int v_idx) {
    switch (v_idx) {
        case 0: return {Side::B, Side::L};
        case 1: return {Side::R, Side::B};
        case 2: return {Side::L, Side::T};
        case 3: return {Side::T, Side::R};
        default: throw std::runtime_error("Invalid base node index");
    }
}

BlockStructured::BlockStructured() : Generator("block-structured") {
    m_fixed = [](const Vector3d&) -> bool {
        return false;
    };
}

Block &BlockStructured::operator[](int idx) const {
    if (m_stage != EDITABLE) {
        throw std::runtime_error("BlockStructured::operator[] is available only at edit stage");
    }
    if (idx >= size() || idx < -size()) {
        throw std::out_of_range("BlockStructured::operator[]: index > blocks.size()");
    }
    if (idx < 0) {
        idx += size();
    }
    if (!m_blocks[idx]) {
        throw std::runtime_error("BlockStructured::operator[]: block does not exist");
    }
    return *m_blocks[idx];
}

void BlockStructured::operator+=(const std::array<BaseNode::Ptr, 4>& base_nodes) {
    Block::Ptr block = Block::create(base_nodes);
    block->set_index(size());
    m_blocks.emplace_back(block);

    // Добавить рёбра
    for (Side side: sides_2D) {
        auto v1 = block->base_node(side, 0);
        auto v2 = block->base_node(side, 1);
        BaseEdge edge{v1, v2};

        if (m_edges.contains(edge)) {
            if (m_edges[edge].b1.expired()) {
                m_edges[edge].b1 = block;
                m_edges[edge].side1 = side;
            }
            else {
                if (m_edges[edge].b1.lock() == block) {
                    std::cerr << "Attempt to add block second time\n";
                }
                else {
                    if (m_edges[edge].b2.expired()) {
                        m_edges[edge].b2 = block;
                        m_edges[edge].side2 = side;
                    }
                    else {
                        if (m_edges[edge].b2.lock() == block) {
                            std::cerr << "Attempt to add block second time\n";
                        }
                        else {
                            std::cerr << "Edge with more than two blocks\n";
                        }
                    }
                }
            }
        }
        else {
            m_edges[edge] = block_pair_t{.b1=block, .side1=side};
        }
    }
}

void BlockStructured::set_boundary(BaseNode::Ref v1, BaseNode::Ref v2, Curve::Ref curve) {
    BaseEdge edge{v1, v2};
    if (m_edges.contains(edge)) {
        if (m_edges[edge].boundary()) {
            if (m_edges[edge].b1.expired()) {
                std::cerr << "No blocks on edge\n";
            }
            else {
                m_edges[edge].b1.lock()->set_boundary(m_edges[edge].side1, curve);
            }
        }
        else {
            std::cerr << "Attempt to set boundary for inner edge\n";
        }
    }
    else {
        throw std::runtime_error("Has no edge for pair of base nodes");
    }
}

void BlockStructured::set_boundary(std::initializer_list<BaseNode::Ptr> nodes, Curve::Ref curve) {
    if (nodes.size() < 2) {
        throw std::invalid_argument("At least two nodes required");
    }

    auto it = nodes.begin();
    BaseNode::Ptr prev = *it++;

    while (it != nodes.end()) {
        BaseNode::Ref current = *it++;
        set_boundary(prev, current, curve);
        prev = current;
    }
}

void BlockStructured::remove_redundant() {
    const auto it = std::ranges::remove(m_blocks, nullptr).begin();
    m_blocks.erase(it, m_blocks.end());
    std::erase_if(m_edges, [](const auto& pair) {
        return pair.second.b1.expired() || pair.second.b2.expired();
    });
    for (int i = 0; i < size(); ++i) {
        m_blocks[i]->set_index(i);
    }
}

void BlockStructured::link_blocks() {
    if (m_stage != EDITABLE) {
        throw std::runtime_error("BlockStructured::link() is available only at edit stage");
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

    // TODO: Проверить связность области!

    m_stage = LINKED;
}

inline std::tuple<double, double> figsize(Box bbox) {
    constexpr double max_width{10.0};
    constexpr double max_height{6.0};

    double aspect = bbox.size()[0] / bbox.size()[1];

    double width{max_width}, height{max_height};
    if (aspect > 1.0) {
        height = width / aspect;
    }
    else {
        width = height * aspect;
    }
    if (width > max_width) {
        height *= max_width / width;
        width = max_width;
    }
    if (height > max_height) {
        width *= max_height / height;
        height = max_height;
    }
    return {width, height};
}

void BlockStructured::plot() const {
    // bounding box
    Box bbox = Box::Empty(2);
    for (const auto& block: m_blocks) {
        auto& vertices = block->mapping();
        for (int i = 0; i < vertices.size1(); ++i) {
            for (int j = 0; j < vertices.size2(); ++j) {
                bbox.capture(vertices(i, j));
            }
        }
    }

    utils::pyplot plt;
    plt.figure({.figsize = figsize(bbox), .dpi = 175});
    plt.set_aspect_equal();

    std::string line_color = "tab:blue";
    std::string outline_color = "black";

    // Сетка внутри блоков
    std::vector<double> xs, ys;
    for (const auto& block: m_blocks) {
        auto& vertices = block->mapping();
        if (vertices.empty()) { continue; }

        xs.resize(vertices.size1());
        ys.resize(vertices.size1());
        for (int j = 1; j < vertices.size2() - 1; ++j) {
            for (int i = 0; i < vertices.size1(); ++i) {
                xs[i] = vertices(i, j).x();
                ys[i] = vertices(i, j).y();
            }
            plt.plot(xs, ys, {.linewidth=1.0, .color=line_color});
        }

        xs.resize(vertices.size2());
        ys.resize(vertices.size2());
        for (int i = 1; i < vertices.size1() - 1; ++i) {
            for (int j = 0; j < vertices.size2(); ++j) {
                xs[j] = vertices(i, j).x();
                ys[j] = vertices(i, j).y();
            }
            plt.plot(xs, ys, {.linewidth=1.0, .color=line_color});
        }
    }

    // Границы блоков
    for (const auto& block: m_blocks) {
        auto& vertices = block->mapping();
        if (vertices.empty()) { continue; }

        xs.resize(vertices.size1());
        ys.resize(vertices.size1());
        for (int j: {0, vertices.size2() - 1}) {
            for (int i = 0; i < vertices.size1(); ++i) {
                xs[i] = vertices(i, j).x();
                ys[i] = vertices(i, j).y();
            }
            plt.plot(xs, ys, {.linewidth=1.5, .color=outline_color});
        }

        xs.resize(vertices.size2());
        ys.resize(vertices.size2());
        for (int i: {0, vertices.size1() - 1}) {
            for (int j = 0; j < vertices.size2(); ++j) {
                xs[j] = vertices(i, j).x();
                ys[j] = vertices(i, j).y();
            }
            plt.plot(xs, ys, {.linewidth=1.5, .color=outline_color});
        }

        // Угловые вершины
        xs.resize(4);
        ys.resize(4);
        for (int i = 0; i < 4; ++i) {
            xs[i] = vertices.corner(i).x();
            ys[i] = vertices.corner(i).y();
        }
        plt.plot(xs, ys, {.linestyle="none", .color=outline_color, .marker="o", .markersize=4});
    }

    plt.tight_layout();
    plt.show();
}

void BlockStructured::plot(const Tables2D& all_vertices) {
    // bounding box
    Box bbox = Box::Empty(2);
    for (const auto& vertices: all_vertices) {
        for (int i = 0; i < vertices.size1(); ++i) {
            for (int j = 0; j < vertices.size2(); ++j) {
                bbox.capture(vertices(i, j)->pos);
            }
        }
    }

    utils::pyplot plt;
    plt.figure({.figsize = figsize(bbox), .dpi = 175});
    plt.set_aspect_equal();

    std::string line_color = "tab:blue";
    std::string outline_color = "black";

    // Сетка внутри блоков
    std::vector<double> xs, ys;
    for (const auto& vertices: all_vertices) {
        if (vertices.empty()) { continue; }

        xs.resize(vertices.size1());
        ys.resize(vertices.size1());
        for (int j = 1; j < vertices.size2() - 1; ++j) {
            for (int i = 0; i < vertices.size1(); ++i) {
                xs[i] = vertices(i, j)->x();
                ys[i] = vertices(i, j)->y();
            }
            plt.plot(xs, ys, {.linewidth=1.0, .color=line_color});
        }

        xs.resize(vertices.size2());
        ys.resize(vertices.size2());
        for (int i = 1; i < vertices.size1() - 1; ++i) {
            for (int j = 0; j < vertices.size2(); ++j) {
                xs[j] = vertices(i, j)->x();
                ys[j] = vertices(i, j)->y();
            }
            plt.plot(xs, ys, {.linewidth=1.0, .color=line_color});
        }
    }

    // Границы блоков
    for (const auto& vertices: all_vertices) {
        if (vertices.empty()) { continue; }

        xs.resize(vertices.size1());
        ys.resize(vertices.size1());
        for (int j: {0, vertices.size2() - 1}) {
            for (int i = 0; i < vertices.size1(); ++i) {
                xs[i] = vertices(i, j)->x();
                ys[i] = vertices(i, j)->y();
            }
            plt.plot(xs, ys, {.linewidth=1.5, .color=outline_color});
        }

        xs.resize(vertices.size2());
        ys.resize(vertices.size2());
        for (int i: {0, vertices.size1() - 1}) {
            for (int j = 0; j < vertices.size2(); ++j) {
                xs[j] = vertices(i, j)->x();
                ys[j] = vertices(i, j)->y();
            }
            plt.plot(xs, ys, {.linewidth=1.5, .color=outline_color});
        }

        // Угловые вершины
        xs.resize(4);
        ys.resize(4);
        for (int i = 0; i < 4; ++i) {
            xs[i] = vertices.corner(i)->x();
            ys[i] = vertices.corner(i)->y();
        }
        plt.plot(xs, ys, {.linestyle="none", .color=outline_color, .marker="o", .markersize=4});
    }

    plt.tight_layout();
    plt.show();
}

void BlockStructured::plot_debug() const {
    utils::pyplot plt;
    plt.figure({.dpi=210});

    plt.set_aspect_equal();

    // Собираем уникальные вершины
    std::set<BaseNode::Ptr> unique_nodes;
    for (const auto& block: m_blocks) {
        for (const auto& node: block->base_nodes()) {
            unique_nodes.insert(node);
        }
    }

    // Построить узлы
    std::vector<double> nodes_x, nodes_y;
    for (const auto& node: unique_nodes) {
        nodes_x.push_back(node->x());
        nodes_y.push_back(node->y());
    }
    //plt.plot(nodes_x, nodes_y, {.linestyle="none", .marker="o"});

    // Проставить смежность
    for (const auto& node: unique_nodes) {
        int n_nodes = node->n_adjacent_nodes();
        int n_blocks = node->n_adjacent_blocks();
        //plt.text(node->x(), node->y(), std::format("{}, {}", n_blocks, n_nodes));
    }

    // Обходы по смежным
    for (const auto& node: unique_nodes) {
        double rn{0.03}, rb{0.05};
        std::vector<double> points_x, points_y;
        for (auto& adj_node: node->adjacent_nodes()) {
            double phi = std::atan2(adj_node.lock()->y() - node->y(), adj_node.lock()->x() - node->x());
            points_x.push_back(node->x() + rn * std::cos(phi));
            points_y.push_back(node->y() + rn * std::sin(phi));
        }
        //plt.plot({points_x.front()}, {points_y.front()}, {.linewidth=3.0, .color="green", .marker="o"});
        //plt.plot(points_x, points_y, {.linewidth=1.0, .color="green", .marker="."});
        //plt.plot({points_x.back()}, {points_y.back()}, {.linewidth=3.0, .color="green", .marker="x"});

        points_x.clear(); points_y.clear();
        for (auto& adj_block: node->adjacent_blocks()) {
            Vector3d vc = adj_block.lock()->center() - node->pos();
            double phi = std::atan2(vc.y(), vc.x());
            points_x.push_back(node->x() + rb * std::cos(phi));
            points_y.push_back(node->y() + rb * std::sin(phi));
        }
        //plt.plot({points_x.front()}, {points_y.front()}, {.linewidth=1.0, .color="blue", .marker="o"});
        //plt.plot(points_x, points_y, {.linewidth=1.0, .color="blue", .marker="."});
        //plt.plot({points_x.back()}, {points_y.back()}, {.linewidth=1.0, .color="blue", .marker="x"});
    }

    // Очертания блоков
    for (const auto& block: m_blocks) {
        std::vector<double> bxs, bys;
        for (int i: {0, 1, 3, 2, 0}) {
            bxs.push_back(block->base_node(i)->x());
            bys.push_back(block->base_node(i)->y());
        }
        //plt.plot(bxs, bys, {.linestyle="solid", .marker="."});
    }

    // Центры блоков
    for (int ib = 0; ib < m_blocks.size(); ++ib) {
        Vector3d bc = m_blocks[ib]->center();
        //plt.text(bc.x(), bc.y(), std::format("B{}", ib));

        //auto ptr1 = reinterpret_cast<uintptr_t>(m_blocks[ib]->lambda_ptr(Axis::X));
        //auto ptr2 = reinterpret_cast<uintptr_t>(m_blocks[ib]->lambda_ptr(Axis::Y));
        //std::cout << std::format("  Block {}: {}, {}\n", ib, ptr1 % 97, ptr2 % 97);
    }

    // Связи блоков
    for (const auto& block: m_blocks) {
        Vector3d c1 = block->center();
        for (Side side: sides_2D) {
            auto adj = block->adjacent_block(side);
            if (!adj) continue;

            Vector3d v1 = block->base_node(side, 0)->pos();
            Vector3d v2 = block->base_node(side, 1)->pos();
            Vector3d c2 = adj->center();

            Vector3d vc = 0.5 * (v1 + v2);
            Vector3d dc = 0.05 * (c2 - c1);

            //plt.arrow(vc.x(), vc.y(), dc.x(), dc.y());
        }
    }

    // Сетка для оптимизации
    for (const auto& block: m_blocks) {
        auto& verts = block->mapping();

        if (verts.empty()) { continue; }
        std::vector<double> xs1(verts.size1());
        std::vector<double> ys1(verts.size1());
        for (int j = 0; j < verts.size2(); ++j) {
            for (int i = 0; i < verts.size1(); ++i) {
                xs1[i] = verts(i, j).x();
                ys1[i] = verts(i, j).y();
            }
            plt.plot(xs1, ys1, {.linewidth=1.0, .color="black"});
        }

        std::vector<double> xs2(verts.size2());
        std::vector<double> ys2(verts.size2());
        for (int i = 0; i < verts.size1(); ++i) {
            for (int j = 0; j < verts.size2(); ++j) {
                xs2[j] = verts(i, j).x();
                ys2[j] = verts(i, j).y();
            }
            plt.plot(xs2, ys2, {.linewidth=1.0, .color="black"});
        }

        //for (int i = 0; i <= size(block, Axis::X); ++i) {
        //    for (int j = 0; j <= size(block, Axis::Y); ++j) {
                // plt.text(verts(i, j).x(), verts(i, j).y(), std::to_string(verts(i,j)->degree()));
        //    }
        //}

        /*
        // Связи с соседями
        for (int i = 0; i <= size(block, Axis::X); ++i) {
            for (int j = 0; j <= size(block, Axis::Y); ++j) {
                Vector3d v1 = verts(i, j)->pos;
                for (auto edge: verts(i, j)->adjacent()) {
                    Vector3d v2 = edge.pos();
                    Vector3d dr = 0.25 * (v2 - v1);

                    //plt.arrow(v1.x(), v1.y(), dr.x(), dr.y());

                    //Vector3d vc = 0.5 * (v1 + v2);
                    //auto ptr1 = reinterpret_cast<uintptr_t>(edge.lambda1);
                    //auto ptr2 = reinterpret_cast<uintptr_t>(edge.lambda2);
                    //if (ptr1 > ptr2) std::swap(ptr1, ptr2);
                    //plt.text(vc.x(), vc.y(), std::format("({}, {})", ptr1 % 97, ptr2 % 97), {.ha = "center", .va = "center"});
                }
            }
        }
        */
    }

    plt.tight_layout();
    plt.show();
}

int BlockStructured::calc_cells(const Pairs<int>& sizes) const {
    int count = 0;
    for (int b1 = 0; b1 < size(); ++b1) {
        count += sizes[b1][Axis::X] * sizes[b1][Axis::Y];
    }
    return count;
}

int BlockStructured::calc_nodes(const Pairs<int>& sizes) const {
    int count = 0;
    for (int b1 = 0; b1 < size(); ++b1) {
        count += (sizes[b1][Axis::X] + 1) * (sizes[b1][Axis::Y] + 1);
    }
    return count;
}

Pairs<int> BlockStructured::auto_block_sizes(int N) {
    for (const auto& block: m_blocks) {
        if (bad_number(block->modulus())) {
            throw std::runtime_error("Need to define modulus before");
        }
    }

    // Сбросить относительные размеры
    Pairs<double> rel_sizes(m_blocks.size());
    for (int b1 = 0; b1 < size(); ++b1) {
        rel_sizes[b1][Axis::X] = NAN;
        rel_sizes[b1][Axis::Y] = NAN;
    }

    // Задаем у первого блока
    init_rel_sizes(rel_sizes, 0);

    // Алгоритм не сработает для несвязных областей
    // TODO: Проверка связности области
    std::list<Block::Ptr> list;
    for (int i = 1; i < m_blocks.size(); ++i) {
        list.emplace_back(m_blocks[i]);
    }

    int count = 0;
    while (!list.empty() && count < 10 * m_blocks.size()) {
        const auto& block = list.front();
        list.pop_front();

        // Размеры блока не определены полностью, переместить в начало
        bool defined = update_rel_sizes(rel_sizes, block);
        if (!defined) {
            list.push_back(block);
        }
        ++count;
    }
    if (count > 10 * m_blocks.size()) {
        throw std::runtime_error("Возможно, не связная область, бесконечный цикл");
    }

    // Находим минимальный относительный размер
    double min_size = std::numeric_limits<double>::max();
    for (int b1 = 0; b1 < size(); ++b1) {
        min_size = std::min(min_size, rel_sizes[b1][Axis::X]);
        min_size = std::min(min_size, rel_sizes[b1][Axis::Y]);
    }

    // Проставить реальные размеры (N ячеек на min_size)
    Pairs<int> sizes(m_blocks.size());
    for (int b1 = 0; b1 < size(); ++b1) {
        sizes[b1][Axis::X] = static_cast<int>(std::round(N * rel_sizes[b1][Axis::X] / min_size));
        sizes[b1][Axis::Y] = static_cast<int>(std::round(N * rel_sizes[b1][Axis::Y] / min_size));
    }
    return sizes;
}

void BlockStructured::conformal_info(const Pairs<double>& lambda) const {
    for (int b1 = 0; b1 <= m_blocks.size(); ++b1) {
        std::cout << std::format("      Block {:3}.  K: {:.2f};  λ_1: {:.2f}; λ_2: {:.2f}",
            b1, m_blocks[b1]->modulus(), lambda[b1][Axis::X], lambda[b1][Axis::Y]);
    }
}

void BlockStructured::optimize(const optimize_options& opts) {
    int n_steps = std::max(1, std::min(opts.steps, 7));
    int N = std::max(2, opts.N);

    // На последней итерации N не должно превышать 1024
    N = std::min(N, static_cast<int>(std::floor(std::pow(2, 11 - n_steps))));

    int verbose = std::max(0, opts.verbose);
    double eps = std::max(1.0e-6, std::min(opts.eps, 0.01));

    if (verbose > 0) {
        std::cout << "Start BlockStructured optimization ";
        std::cout << std::format("({} steps; start with {} cells; eps: {:.2e}; verbose: {})\n", n_steps, N, eps, verbose);
    }
    for (int step = 0; step < n_steps; ++step) {
        if (verbose > 0) {
            std::cout << std::format("Optimize step {}. N: {}\n", step + 1, N);
        }
        optimize(N, eps, verbose);
        N *= 2;
        eps /= 2.0;
    }
}

void BlockStructured::optimize(int N, double eps, int verbose) {
    if (m_blocks.empty()) {
        throw std::runtime_error("Add blocks before optimization");
    }
    if (m_stage == EDITABLE) {
        if (verbose > 1) {
            std::cout << "  Link blocks\n";
        }
        link_blocks();
        if (verbose > 1) {
            std::cout << "  Init optimization\n";
        }
        optimize_init(N, eps, verbose);
        return;
    }
    if (m_stage == OPTIMIZED) {
        if (verbose > 1) {
            std::cout << "  Repeat optimization\n";
        }
        optimize_again(N, eps, verbose);
        return;
    }
}

inline AxisPair<double> calc_lambda_zero(double K, AxisPair<int> sizes) {
    AxisPair<double> lambda{0.0, 0.0};
    if (sizes[Axis::X] >= 1 && sizes[Axis::Y] >= 1) {
        lambda[Axis::Y] = K * sizes[Axis::Y] / sizes[Axis::X];
        lambda[Axis::X] = 1.0 / lambda[Axis::Y];
    }
    return lambda;
}

inline AxisPair<double> calc_lambda_nan(double K, AxisPair<int> sizes) {
    AxisPair<double> lambda{NAN, NAN};
    if (good_number(K) && sizes[Axis::X] >= 1 && sizes[Axis::Y] >= 1) {
        lambda[Axis::Y] = K * sizes[Axis::Y] / sizes[Axis::X];
        lambda[Axis::X] = 1.0 / lambda[Axis::Y];
    }
    return lambda;
}

Pairs<double> BlockStructured::get_lambda_nan(const Pairs<int>& sizes) const {
    if (m_blocks.size() != sizes.size()) {
        throw std::runtime_error("get_lambda_nan: block size mismatch");
    }
    Pairs<double> lambda(m_blocks.size());
    for (int b1 = 0; b1 < size(); ++b1) {
        lambda[b1] = calc_lambda_nan(m_blocks[b1]->modulus(), sizes[b1]);
    }
    return lambda;
}

void BlockStructured::optimize_init(int N, double eps, int verbose) {
    // Оценить конформный модуль по геометрии
    for (int b1 = 0; b1 < size(); ++b1) {
        double K = m_blocks[b1]->estimate_modulus();
        m_blocks[b1]->set_modulus(K);
    }

    // Проставить размеры блоков
    const auto sizes = auto_block_sizes(N);
    auto lambda = get_lambda_nan(sizes);

    if (verbose > 2) {
        sizes_info(sizes);
    }
    if (verbose > 1) {
        std::cout << "    Cells count: " << calc_cells(sizes) << "\n";
    }

    // Первое создание вершин
    check_consistency(sizes);
    auto vertices = create_vertices(sizes);
    merge_vertices(vertices, sizes);
    link_vertices(vertices, sizes, lambda);
    if (verbose > 4) {
        plot(vertices);
    }

    int counter = 0;
    double error = 1.0;

    while (error > eps && counter < 5000) {
        if (verbose > 2 && counter % 50 == 0) {
            std::cout << std::format("    step {:6}, eps: {:.2e}\n", counter, error);
            if (verbose > 3 && counter % 50 == 0) {
                conformal_info(lambda);
            }
        }
        smooth_vertices(vertices);
        error = update_vertices(vertices);
        for (int b1 = 0; b1 < size(); ++b1) {
            m_blocks[b1]->update_modulus(vertices[b1]);
            lambda[b1] = calc_lambda_zero(m_blocks[b1]->modulus(), sizes[b1]);
        }
        ++counter;
    }
    if (verbose > 2) {
        std::cout << std::format("    step {:6}, eps: {:.2e}\n", counter, error);
        if (verbose > 3) {
            conformal_info(lambda);
        }
    }
    if (verbose > 4) {
        plot(vertices);
    }
    for (int k = 0; k < size(); ++k) {
        m_blocks[k]->set_mapping(vertices[k]);
    }
    m_stage = OPTIMIZED;
}

void BlockStructured::optimize_again(int N, double eps, int verbose) {
    // Проставить размеры блоков
    const auto sizes = auto_block_sizes(N);
    auto lambda = get_lambda_nan(sizes);

    if (verbose > 2) {
        sizes_info(sizes);
    }
    if (verbose > 1) {
        std::cout << "    Cells count: " << calc_cells(sizes) << "\n";
    }

    // Повторное создание вершин
    check_consistency(sizes);

    auto vertices = create_vertices_again(sizes);
    merge_vertices(vertices, sizes);
    link_vertices(vertices, sizes, lambda);

    int counter = 0;
    double error = 1.0;

    while (error > eps && counter < 5000) {
        if (verbose > 2 && counter % 50 == 0) {
            std::cout << std::format("    step {:6}, eps: {:.2e}\n", counter, error);
            if (verbose > 3) {
                conformal_info(lambda);
            }
        }
        smooth_vertices(vertices);
        error = update_vertices(vertices);
        for (int b1 = 0; b1 < size(); ++b1) {
            m_blocks[b1]->update_modulus(vertices[b1]);
            lambda[b1] = calc_lambda_zero(m_blocks[b1]->modulus(), sizes[b1]);
        }
        ++counter;
    }
    if (verbose > 2) {
        std::cout << std::format("    step {:6}, eps: {:.2e}\n", counter, error);
        if (verbose > 3) {
            conformal_info(lambda);
        }
    }
    if (verbose > 4) {
        plot(vertices);
    }
    for (int k = 0; k < size(); ++k) {
        m_blocks[k]->set_mapping(vertices[k]);
    }
    m_stage = OPTIMIZED;
}

Grid BlockStructured::make() {
    throw std::runtime_error("Not implemented");
    Grid grid;
    //grid.reserve_cells(calc_cells());
    //grid.reserve_nodes(calc_nodes());

    /*

    // Ставим неопределенными
    for (auto &block: m_blocks) {
        for (int i = 0; i <= block.size1(); ++i) {
            for (int j = 0; j <= block.size2(); ++j) {
                block.vertex(i, j)->index = -1;
            }
        }
    }

    int n_nodes = 0;
    for (auto &block: m_blocks) {
        for (int i = 0; i <= block.size1(); ++i) {
            for (int j = 0; j <= block.size2(); ++j) {
                BsVertex::Ref v = block.vertex(i, j);
                if (v->index < 0) {
                    v->index = n_nodes;

                    GNode::Ptr node = GNode::create(v->v1);
                    node->index = n_nodes;
                    node->set_boundaries(v->boundaries());
                    grid += node;
                    ++n_nodes;
                }
            }
        }
    }

    int n_cells = 0;
    for (auto &block: m_blocks) {
        for (int i = 0; i < block.size1(); ++i) {
            for (int j = 0; j < block.size2(); ++j) {
                GCell cell = GCell::quad(
                        {
                           grid.node(block.vertex(i, j)->index),
                           grid.node(block.vertex(i+1, j)->index),
                           grid.node(block.vertex(i+1, j+1)->index),
                           grid.node(block.vertex(i, j+1)->index)
                        });
                cell.index = n_cells;
                grid += cell;
                ++n_cells;
            }
        }
    }
    
    grid.setup_adjacency();
    */

    return grid;
}

void BlockStructured::check_consistency(const Pairs<int>& sizes) const {
    for (int b1 = 0; b1 < size(); ++b1) {
        auto block = m_blocks[b1];

        if (sizes[b1][Axis::X] < 1 || sizes[b1][Axis::Y] < 1) {
            throw std::runtime_error("Zero size of some block");
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
                    std::string message = std::format("Block::check_consistency error: size mismatch (blocks {}, {})", b1, b2);
                    std::cerr << message << "\n";
                    throw std::runtime_error(message);
                }
            }
            else {
                std::string message = "Block::check_consistency error: Boundary or neighbor should be defined";
                std::cerr << message << "\n";
                throw std::runtime_error(message);
            }
        }
    }
}

Tables2D BlockStructured::create_vertices(const Pairs<int>& sizes) const {
    if (m_blocks.size() != sizes.size()) {
        throw std::runtime_error("BlockStructured::create_vertices: block size mismatch");
    }
    Tables2D vertices(m_blocks.size());
    for (int b1 = 0; b1 < size(); ++b1) {
        vertices[b1] = m_blocks[b1]->create_vertices(sizes[b1]);
    }
    return vertices;
}

Tables2D BlockStructured::create_vertices_again(const Pairs<int>& sizes) const {
    if (m_blocks.size() != sizes.size()) {
        throw std::runtime_error("BlockStructured::create_vertices: block size mismatch");
    }
    Tables2D vertices(m_blocks.size());
    for (int b1 = 0; b1 < size(); ++b1) {
        vertices[b1] = m_blocks[b1]->create_vertices_again(sizes[b1]);
    }
    return vertices;
}

void BlockStructured::merge_vertices(Tables2D& all_vertices, const Pairs<int>& sizes) const {
    for (int b1 = 0; b1 < size(); ++b1) {
        auto block = m_blocks[b1];

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
}

void BlockStructured::link_vertices(Tables2D& all_vertices, const Pairs<int>& sizes, Pairs<double>& lambda) const {
    for (int b1 = 0; b1 < size(); ++b1) {
        auto self = m_blocks[b1];
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
                throw std::runtime_error("BFace is not boundary and has no neighbor");
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
                    throw std::runtime_error("Face is not boundary and has no neighbor");
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
            if (node->is_fixed()) { continue; }

            // Угловая вершина области, должна быть неподвижной
            if (node->n_adjacent_blocks() < 2) { continue; }

            // Угловая вершина области по другому критерию
            if (auto [side1, side2] = node_adjacent_sides(v_idx);
                self->boundary(side1) && self->boundary(side2)) {
                continue;
            }

            // Вершина на границе
            if (node->is_boundary()) {
                // Предполагаем, что смежные вершины и блоки упорядочены верно,
                // то есть против часовой стрелки внутри области
                auto& adjacent_nodes = node->adjacent_nodes();
                auto& adjacent_blocks = node->adjacent_blocks();
                if (adjacent_blocks.size() + 1 != adjacent_nodes.size()) {
                    throw std::runtime_error("Block::link_vertices: n_blocks + 1 != n_nodes for boundary node");
                }

                auto corner = vertices.corner(v_idx);
                std::vector<BsEdge> edges;

                // Первая граничная вершина
                {
                    BaseNode::Ptr node2 = adjacent_nodes.front().lock();
                    Block::Ptr block = adjacent_blocks.front().lock();
                    if (!node2 || !block) { throw std::runtime_error("Block::link_vertices: invalid block"); }
                    auto side = block->get_side(node, node2);
                    if (!block->boundary(side)) {
                        throw std::runtime_error("Block::link_vertices: invalid block (not boundary)");
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
                    if (!node2) { throw std::runtime_error("Block::link_vertices: invalid node"); }
                    Block::Ptr block = adjacent_blocks[i].lock();
                    auto prev_block = adjacent_blocks[i - 1].lock();
                    if (!block || !prev_block) { throw std::runtime_error("Block::link_vertices: invalid block"); }
                    int b3 = block->index();
                    Side side = block->get_side(node, node2);
                    Side prev_side = prev_block->get_side(node, node2);
                    edges.push_back(BsEdge::Inside(
                        all_vertices[b3].boundary(side, 1),
                        &lambda[b3][side],
                        &lambda[prev_block->index()][prev_side]
                    ));
                    if (block->boundary(side)) {
                        throw std::runtime_error("Block::link_vertices: invalid block (why boundary)");
                    }
                }
                // Последняя граничная вершина
                {
                    BaseNode::Ptr node2 = adjacent_nodes.back().lock();
                    Block::Ptr block = adjacent_blocks.back().lock();
                    if (!node2 || !block) { throw std::runtime_error("Block::link_vertices: invalid block"); }
                    Side side = block->get_side(node, node2);
                    edges.push_back(BsEdge::Border(
                        all_vertices[block->index()].boundary(side, sizes[block->index()][side] - 1),
                        &lambda[block->index()][side]
                    ));
                    if (!block->boundary(side)) {
                        throw std::runtime_error("Block::link_vertices: invalid block (not boundary)");
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
                    throw std::runtime_error("Block::link_vertices: n_blocks != n_nodes for inner node");
                }

                auto corner = vertices.corner(v_idx);
                std::vector<BsEdge> edges;

                for (int i = 0; i < adjacent_blocks.size(); ++i) {
                    BaseNode::Ptr node2 = adjacent_nodes[i].lock();
                    if (!node2) { throw std::runtime_error("Block::link_vertices: invalid node"); }
                    Block::Ptr block = adjacent_blocks[i].lock();
                    int j = (i + adjacent_blocks.size() - 1) % adjacent_blocks.size();
                    auto prev_block = adjacent_blocks[j].lock();
                    if (!block || !prev_block) { throw std::runtime_error("Block::link_vertices: invalid block"); }
                    Side side = block->get_side(node, node2);
                    Side prev_side = prev_block->get_side(node, node2);
                    edges.push_back(BsEdge::Inside(
                        all_vertices[block->index()].boundary(side, 1),
                        &lambda[block->index()][side],
                        &lambda[prev_block->index()][prev_side]
                    ));
                    if (block->boundary(side)) {
                        throw std::runtime_error("Block::link_vertices: invalid block (why boundary)");
                    }
                }
                corner->set_edges(edges);
            }
        }
    }
}

inline void smooth_one(const Array2D<BsVertex::Ptr>& vertices) {
    for (int i = 0; i < vertices.size(Axis::X); ++i) {
        for (int j = 0; j < vertices.size(Axis::Y); ++j) {
            auto vc = vertices(i, j);
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

void BlockStructured::smooth_vertices(const Tables2D& all_vertices) {
    for (const auto& vertices: all_vertices) {
        smooth_one(vertices);
    }
}

double BlockStructured::update_vertices(const Tables2D& all_vertices) {
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

// ---------- Размеры блока ----------------------------------------- --------------------------------------------------

// Проставить размеры, в том числе у связанных блоков
template <typename Type>
void SET_SIZE(const std::vector<Block::Ptr>& blocks, Pairs<Type>& sizes, int b1, Axis axis, Type N) {
    if (blocks.size() != sizes.size()) {
        throw std::runtime_error("SET SIZE: SIZE MISMATCH");
    }
    if (sizes[b1][axis] == N) { return; }
    sizes[b1][axis] = N;
    for (Side side: sides_by_axis(axis)) {
        auto neib = blocks[b1]->adjacent_block(side);
        if (neib) {
            int b2 = neib->index();
            auto twin = blocks[b1]->twin_face(side);
            if (sizes[b2][twin] != N) {
                SET_SIZE(blocks, sizes, b2, to_axis(twin), N);
            }
        }
    }
}

void BlockStructured::set_size(Pairs<int>& sizes, int b1, Axis axis, int N) const {
    SET_SIZE(m_blocks, sizes, b1, axis, N);
}

void BlockStructured::set_rel_size(Pairs<double>& sizes, int b1, Axis axis, double N) const {
    SET_SIZE(m_blocks, sizes, b1, axis, N);
}

void BlockStructured::init_rel_sizes(Pairs<double>& sizes, int b1) const {
    if (std::isnan(m_blocks[b1]->modulus())) {
        throw std::runtime_error("Invalid Modulus");
    }
    set_rel_size(sizes, b1, Axis::X, m_blocks[b1]->modulus());
    set_rel_size(sizes, b1, Axis::Y, 1.0);
}

bool BlockStructured::update_rel_sizes(Pairs<double>& sizes, Block::Ref b) const {
    int b1 = b->index();
    if (std::isnan(b->modulus())) {
        throw std::runtime_error("update_rel_sizes error: Invalid Modulus");
    }
    if (bad_number(sizes[b1][Axis::X]) && good_number(sizes[b1][Axis::Y])) {
        set_rel_size(sizes, b1, Axis::X, sizes[b1][Axis::Y] * b->modulus());
    }
    if (bad_number(sizes[b1][Axis::Y]) && good_number(sizes[b1][Axis::X])) {
        set_rel_size(sizes, b1, Axis::Y, sizes[b1][Axis::X] / b->modulus());
    }
    return good_number(sizes[b1][Axis::X]) && good_number(sizes[b1][Axis::Y]);
}

void BlockStructured::sizes_info(const Pairs<int>& sizes) const {
    for (int b1 = 0; b1 < size(); ++b1) {
        std::cout << std::format("    Block {:3}.  K: {:.2f};  sizes: ({:3}, {:3})", b1,
            m_blocks[b1]->modulus(), sizes[b1][Axis::X], sizes[b1][Axis::Y]);
    }
}

} // namespace zephyr::geom::generator
