#include <array>
#include <algorithm>
#include <format>
#include <map>
#include <iomanip>
#include <list>
#include <ranges>
#include <bits/fs_fwd.h>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>
#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/block_structured.h>
#include <zephyr/utils/pyplot.h>

namespace zephyr::geom::generator {

constexpr int node_idx(Side side, int idx) {
    switch (side) {
        case Side::L: return std::array{0, 2}[idx];
        case Side::R: return std::array{1, 3}[idx];
        case Side::B: return std::array{0, 1}[idx];
        case Side::T: return std::array{2, 3}[idx];
        default: throw std::runtime_error("Invalid side");
    }
}

inline bool bad_number(double num) {
    return std::isinf(num) || std::isnan(num) || num <= 0.0;
}

inline bool good_number(double num) {
    return !std::isinf(num) && !std::isnan(num) && num > 0.0;
}

// Перпендикулярная ось
constexpr Axis opposite(Axis axis) {
    return axis == Axis::X ? Axis::Y : Axis::X;
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
        auto[v1, v2] = block->base_nodes(side);
        node_pair_t edge{v1, v2};

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
    node_pair_t edge{v1, v2};
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
        auto& vertices = block->m_mapping;
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
        auto& vertices = block->m_mapping;
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
        auto& vertices = block->m_mapping;
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

        auto ptr1 = reinterpret_cast<uintptr_t>(m_blocks[ib]->lambda_ptr(Axis::X));
        auto ptr2 = reinterpret_cast<uintptr_t>(m_blocks[ib]->lambda_ptr(Axis::Y));
        //std::cout << std::format("  Block {}: {}, {}\n", ib, ptr1 % 97, ptr2 % 97);
    }

    // Связи блоков
    for (const auto& block: m_blocks) {
        Vector3d c1 = block->center();
        for (Side side: sides_2D) {
            auto adj = block->adjacent_block(side);
            if (!adj) continue;

            Vector3d v1 = block->base_node(node_idx(side, 0))->pos();
            Vector3d v2 = block->base_node(node_idx(side, 1))->pos();
            Vector3d c2 = adj->center();

            Vector3d vc = 0.5 * (v1 + v2);
            Vector3d dc = 0.05 * (c2 - c1);

            //plt.arrow(vc.x(), vc.y(), dc.x(), dc.y());
        }
    }

    // Сетка для оптимизации
    for (const auto& block: m_blocks) {
        auto& verts = block->m_mapping;

        if (verts.empty()) { continue; }
        std::vector<double> xs1(block->size1() + 1);
        std::vector<double> ys1(block->size1() + 1);
        for (int j = 0; j <= block->size2(); ++j) {
            for (int i = 0; i <= block->size1(); ++i) {
                xs1[i] = verts(i, j).x();
                ys1[i] = verts(i, j).y();
            }
            plt.plot(xs1, ys1, {.linewidth=1.0, .color="black"});
        }

        std::vector<double> xs2(block->size2() + 1);
        std::vector<double> ys2(block->size2() + 1);
        for (int i = 0; i <= block->size1(); ++i) {
            for (int j = 0; j <= block->size2(); ++j) {
                xs2[j] = verts(i, j).x();
                ys2[j] = verts(i, j).y();
            }
            plt.plot(xs2, ys2, {.linewidth=1.0, .color="black"});
        }

        for (int i = 0; i <= block->size1(); ++i) {
            for (int j = 0; j <= block->size2(); ++j) {
                // plt.text(verts(i, j).x(), verts(i, j).y(), std::to_string(verts(i,j)->degree()));
            }
        }

        /*
        // Связи с соседями
        for (int i = 0; i <= block->size1(); ++i) {
            for (int j = 0; j <= block->size2(); ++j) {
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

int BlockStructured::calc_cells() const {
    int count = 0;
    for (auto& block: m_blocks) {
        count += block->size1() * block->size2();
    }
    return count;
}

int BlockStructured::calc_nodes() const {
    int count = 0;
    for (auto& block: m_blocks) {
        count += (block->size1() + 1) * (block->size2() + 1);
    }
    return count;
}

void BlockStructured::setup_block_sizes(int N) {
    for (const auto& block: m_blocks) {
        if (bad_number(block->modulus())) {
            throw std::runtime_error("Need to define modulus before");
        }
    }

    // Сбросить относительные размеры
    for (const auto& block: m_blocks) {
        block->reset_rel_sizes();
    }

    // Задаем у первого блока
    m_blocks.front()->init_rel_sizes();

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
        bool defined = block->update_rel_sizes();
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
    for (const auto &block: m_blocks) {
        min_size = std::min(min_size, block->rel_size(Axis::X));
        min_size = std::min(min_size, block->rel_size(Axis::Y));
    }

    // Проставить реальные размеры (N ячеек на min_size)
    for (const auto &block: m_blocks) {
        block->set_size(Axis::X, static_cast<int>(std::round(N * block->rel_size(Axis::X) / min_size)));
        block->set_size(Axis::Y, static_cast<int>(std::round(N * block->rel_size(Axis::Y) / min_size)));
    }

    /*
    for (const auto& block: m_blocks) {
        std::cout << std::format("Block {}; modulus: {:.2f};\trel_sizes: ({:.2f}, {:.2f})",
            block->index(), block->modulus(), block->rel_size(Axis::X), block->rel_size(Axis::Y));
    }
    */
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

void BlockStructured::optimize_init(int N, double eps, int verbose) {
    // Оценить конформный модуль по геометрии
    for (const auto &block: m_blocks) {
        block->estimate_modulus();
    }

    // Проставить размеры блоков
    setup_block_sizes(N);

    if (verbose > 2) {
        for (const auto& block: m_blocks) {
            std::cout << "    " << block->sizes_info() << "\n";
        }
    }
    if (verbose > 1) {
        std::cout << "    Cells count: " << calc_cells() << "\n";
    }

    // Первое создание вершин
    check_consistency();
    std::vector<axis_pair<int>> all_sizes(m_blocks.size());
    for (int k = 0; k < size(); ++k) {
        all_sizes[k][Axis::X] = m_blocks[k]->size(Axis::X);
        all_sizes[k][Axis::Y] = m_blocks[k]->size(Axis::Y);
    }
    auto vertices = create_vertices(all_sizes);
    merge_vertices(vertices);
    link_vertices(vertices);
    if (verbose > 4) {
        plot(vertices);
    }

    int counter = 0;
    double error = 1.0;

    while (error > eps && counter < 5000) {
        if (verbose > 2 && counter % 50 == 0) {
            std::cout << std::format("    step {:6}, eps: {:.2e}\n", counter, error);
            if (verbose > 3 && counter % 50 == 0) {
                for (const auto &block: m_blocks) {
                    std::cout << "      " << block->conformal_info() << "\n";
                }
            }
        }
        smooth_vertices(vertices);
        error = update_vertices(vertices);
        for (int b1 = 0; b1 < size(); ++b1) {
            m_blocks[b1]->update_modulus(vertices[b1]);
        }
        ++counter;
    }
    if (verbose > 2) {
        std::cout << std::format("    step {:6}, eps: {:.2e}\n", counter, error);
        if (verbose > 3) {
            for (const auto &block: m_blocks) {
                std::cout << "      " << block->conformal_info() << "\n";
            }
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
    setup_block_sizes(N);

    if (verbose > 2) {
        for (const auto& block: m_blocks) {
            std::cout << "    " << block->sizes_info() << "\n";
        }
    }
    if (verbose > 1) {
        std::cout << "    Cells count: " << calc_cells() << "\n";
    }

    // Повторное создание вершин
    check_consistency();
    std::vector<axis_pair<int>> all_sizes(m_blocks.size());
    for (int k = 0; k < size(); ++k) {
        all_sizes[k][Axis::X] = m_blocks[k]->size(Axis::X);
        all_sizes[k][Axis::Y] = m_blocks[k]->size(Axis::Y);
    }

    auto vertices = create_vertices_again(all_sizes);
    merge_vertices(vertices);
    link_vertices(vertices);

    int counter = 0;
    double error = 1.0;

    while (error > eps && counter < 5000) {
        if (verbose > 2 && counter % 50 == 0) {
            std::cout << std::format("    step {:6}, eps: {:.2e}\n", counter, error);
            if (verbose > 3) {
                for (const auto &block: m_blocks) {
                    std::cout << "      " << block->conformal_info() << "\n";
                }
            }
        }
        smooth_vertices(vertices);
        error = update_vertices(vertices);
        for (int b1 = 0; b1 < size(); ++b1) {
            m_blocks[b1]->update_modulus(vertices[b1]);
        }
        ++counter;
    }
    if (verbose > 2) {
        std::cout << std::format("    step {:6}, eps: {:.2e}\n", counter, error);
        if (verbose > 3) {
            for (const auto &block: m_blocks) {
                std::cout << "      " << block->conformal_info() << "\n";
            }
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
    grid.reserve_cells(calc_cells());
    grid.reserve_nodes(calc_nodes());

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

void BlockStructured::check_consistency() const {
    for (int b1 = 0; b1 < size(); ++b1) {
        auto block = m_blocks[b1];

        if (block->size1() < 1 || block->size2() < 1) {
            throw std::runtime_error("Zero size of some block");
        }
        for (Side side: sides_2D) {
            if (block->boundary(side)) { continue; }

            auto neib = block->adjacent_block(side);
            if (neib) {
                Side twin = block->twin_face(side);

                int size1 = block->size(side);
                int size2 = neib->size(twin);

                if (size1 != size2) {
                    std::string message = std::format("Block::check_consistency error: size mismatch (blocks {}, {})", b1, neib->index());
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

Tables2D BlockStructured::create_vertices(const std::vector<axis_pair<int>>& sizes) const {
    if (m_blocks.size() != sizes.size()) {
        throw std::runtime_error("BlockStructured::create_vertices: block size mismatch");
    }
    Tables2D vertices(m_blocks.size());
    for (int b1 = 0; b1 < size(); ++b1) {
        vertices[b1] = m_blocks[b1]->create_vertices(sizes[b1]);
    }
    return vertices;
}

Tables2D BlockStructured::create_vertices_again(const std::vector<axis_pair<int>>& sizes) const {
    if (m_blocks.size() != sizes.size()) {
        throw std::runtime_error("BlockStructured::create_vertices: block size mismatch");
    }
    Tables2D vertices(m_blocks.size());
    for (int b1 = 0; b1 < size(); ++b1) {
        vertices[b1] = m_blocks[b1]->create_vertices_again(sizes[b1]);
    }
    return vertices;
}

void BlockStructured::merge_vertices(Tables2D& all_vertices) const {
    for (int b1 = 0; b1 < size(); ++b1) {
        auto block = m_blocks[b1];

        // Связать вершины с соседом
        for (Side side: sides_2D) {
            auto neib = block->adjacent_block(side);

            if (!neib) { continue; }

            // Сосед есть, но без вершин
            if (all_vertices[neib->index()].empty()) {
                continue;
            }

            Side twin = block->twin_face(side);

            int b2 = neib->index();

            int N = block->size(side);
            for (int idx = 0; idx <= N; ++idx) {
                all_vertices[b1].boundary(side, idx) = all_vertices[b2].boundary(twin, N - idx);
            }
        }
    }
}

void BlockStructured::link_vertices(Tables2D& all_vertices) const {
    for (int b1 = 0; b1 < size(); ++b1) {
        auto self = m_blocks[b1];

        // Очистить все списки
        for (int i = 0; i <= self->size1(); ++i) {
            for (int j = 0; j <= self->size2(); ++j) {
                all_vertices[self->index()](i, j)->clear();
            }
        }

        // Внутренние вершины блоков
        for (int i = 1; i < self->size1(); ++i) {
            for (int j = 1; j < self->size2(); ++j) {
                // Справа, сверху, слева, снизу
                std::vector edges {
                    BsEdge::Inside(all_vertices[self->index()](i + 1, j), self->lambda_ptr(Axis::X)),
                    BsEdge::Inside(all_vertices[self->index()](i, j + 1), self->lambda_ptr(Axis::Y)),
                    BsEdge::Inside(all_vertices[self->index()](i - 1, j), self->lambda_ptr(Axis::X)),
                    BsEdge::Inside(all_vertices[self->index()](i, j - 1), self->lambda_ptr(Axis::Y))
                };
                all_vertices[self->index()](i, j)->set_edges(edges);
            }
        }

        // Вершины на границах блока (без угловых)
        for (Side side: sides_2D) {
            if (!self->boundary(side) && !self->adjacent_block(side)) {
                throw std::runtime_error("BFace is not boundary and has no neighbor");
            }

            Axis axis = to_axis(side);
            Axis perp_axis = opposite(axis);

            if (self->boundary(side)) {
                for (int idx = 1; idx < self->size(side); ++idx) {
                    // Добавить границу
                    all_vertices[self->index()].boundary(side, idx)->add_boundary(self->boundary(side).get());
                    // Добавить связи (три штуки)
                    std::vector edges {
                        BsEdge::Border(all_vertices[self->index()].boundary(side, idx + 1), self->lambda_ptr(axis)),
                        BsEdge::Inside(all_vertices[self->index()].near_boundary(side, idx),  self->lambda_ptr(perp_axis)),
                        BsEdge::Border(all_vertices[self->index()].boundary(side, idx - 1), self->lambda_ptr(axis))
                    };
                    all_vertices[self->index()].boundary(side, idx)->set_edges(edges);
                }
            } else {
                auto neib = self->adjacent_block(side);
                if (!neib) {
                    throw std::runtime_error("Face is not boundary and has no neighbor");
                }

                int N = self->size(side);
                Side twin = self->twin_face(side);

                Axis neib_axis = to_axis(twin);
                Axis prep_neib_axis = opposite(neib_axis);
                for (int idx = 1; idx < N; ++idx) {
                    std::vector edges {
                        BsEdge::Inside(all_vertices[self->index()].boundary(side, idx + 1), self->lambda_ptr(axis), neib->lambda_ptr(neib_axis)),
                        BsEdge::Inside(all_vertices[self->index()].near_boundary(side, idx), self->lambda_ptr(perp_axis)),
                        BsEdge::Inside(all_vertices[self->index()].boundary(side, idx - 1), self->lambda_ptr(axis), neib->lambda_ptr(neib_axis)),
                        BsEdge::Inside(all_vertices[neib->index()].near_boundary(twin, N - idx), neib->lambda_ptr(prep_neib_axis))
                    };
                    all_vertices[self->index()].boundary(side, idx)->set_edges(edges);
                }
            }
        }

        // Угловые вершины блока (здесь могут быть сингулярности)
        for (int v_idx = 0; v_idx < 4; ++v_idx) {
            auto node = self->base_node(v_idx);
            all_vertices[self->index()].corner(v_idx)->clear();

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

                auto corner = all_vertices[self->index()].corner(v_idx);
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
                        block->lambda_ptr(side)
                    ));
                }
                // Внутренние вершины
                for (int i = 1; i < adjacent_blocks.size(); ++i) {
                    BaseNode::Ptr node2 = adjacent_nodes[i].lock();
                    if (!node2) { throw std::runtime_error("Block::link_vertices: invalid node"); }
                    Block::Ptr block = adjacent_blocks[i].lock();
                    auto prev_block = adjacent_blocks[i - 1].lock();
                    if (!block || !prev_block) { throw std::runtime_error("Block::link_vertices: invalid block"); }
                    Side side = block->get_side(node, node2);
                    Side prev_side = prev_block->get_side(node, node2);
                    edges.push_back(BsEdge::Inside(
                        all_vertices[block->index()].boundary(side, 1),
                        block->lambda_ptr(side),
                        prev_block->lambda_ptr(prev_side)
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
                        all_vertices[block->index()].boundary(side, block->size(side) - 1),
                        block->lambda_ptr(side)
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

                auto corner = all_vertices[self->index()].corner(v_idx);
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
                        block->lambda_ptr(side),
                        prev_block->lambda_ptr(prev_side)
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

} // namespace zephyr::geom::generator
