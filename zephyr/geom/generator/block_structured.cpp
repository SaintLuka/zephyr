#include <array>
#include <algorithm>
#include <format>
#include <map>
#include <iomanip>
#include <ranges>

#include <zephyr/geom/grid.h>
#include <zephyr/mesh/side.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>
#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/block_structured.h>
#include <zephyr/utils/pyplot.h>


namespace zephyr::geom::generator {

constexpr int node(Side2D side, int idx) {
    switch (side) {
        case Side2D::L: return std::array{0, 2}[idx];
        case Side2D::R: return std::array{1, 3}[idx];
        case Side2D::B: return std::array{0, 1}[idx];
        case Side2D::T: return std::array{2, 3}[idx];
        default: throw std::runtime_error("Invalid side");
    }
}

BlockStructured::BlockStructured(int n_blocks)
    : Generator("block_structured") {

    m_blocks.resize(n_blocks, nullptr);

    m_fixed = [](const Vector3d&) -> bool {
        return false;
    };
}

void BlockStructured::resize(int n_blocks) {
    if (m_stage != EDITABLE) {
        throw std::runtime_error("BlockStructured::resize() is available only at edit stage");
    }
    if (n_blocks < m_blocks.size()) {
        throw std::runtime_error("Can only increase blocks number");
    }
    m_blocks.resize(n_blocks, nullptr);
}

Block &BlockStructured::operator[](int idx) {
    if (m_stage != EDITABLE) {
        throw std::runtime_error("BlockStructured::operator[] is available only at edit stage");
    }
    if (idx >= m_blocks.size()) {
        throw std::out_of_range("BlockStructured::operator[]");
    }
    if (!m_blocks[idx]) { // блок ещё не создан
        m_blocks[idx] = std::make_shared<Block>();
        m_blocks[idx]->index = idx;
    }
    return *m_blocks[idx];
}

void BlockStructured::remove_null_blocks() {
    const auto it = std::remove(m_blocks.begin(), m_blocks.end(), nullptr);
    m_blocks.erase(it, m_blocks.end());
}

void BlockStructured::link() {
    if (m_stage != EDITABLE) {
        throw std::runtime_error("BlockStructured::link() is available only at edit stage");
    }

    // Удалить нулевые блоки
    remove_null_blocks();

    // Собираем уникальные узлы
    std::map<BaseNode::Ptr, std::set<Block::Ptr>> nodes;
    for (const auto& block: m_blocks) {
        for (const auto& v: block->base_nodes()) {
            v->reset();
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

void BlockStructured::plot() const {
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

        auto ptr1 = reinterpret_cast<uintptr_t>(m_blocks[ib]->ratio_ptr(0));
        auto ptr2 = reinterpret_cast<uintptr_t>(m_blocks[ib]->ratio_ptr(1));
        //std::cout << std::format("  Block {}: {}, {}\n", ib, ptr1 % 97, ptr2 % 97);
    }

    // Связи блоков
    for (const auto& block: m_blocks) {
        Vector3d c1 = block->center();
        for (Side2D side: Side2D::items()) {
            auto adj = block->adjacent_block(side);
            if (!adj) continue;

            Vector3d v1 = block->base_node(node(side, 0))->pos();
            Vector3d v2 = block->base_node(node(side, 1))->pos();
            Vector3d c2 = adj->center();

            Vector3d vc = 0.5 * (v1 + v2);
            Vector3d dc = 0.05 * (c2 - c1);

            //plt.arrow(vc.x(), vc.y(), dc.x(), dc.y());
        }
    }

    // Сетка для оптимизации
    for (const auto& block: m_blocks) {
        auto& verts = block->vertices();

        if (verts.empty()) { continue; }
        std::vector<double> xs1(block->size1() + 1);
        std::vector<double> ys1(block->size1() + 1);
        for (int j = 0; j <= block->size2(); ++j) {
            for (int i = 0; i <= block->size1(); ++i) {
                xs1[i] = verts(i, j)->x();
                ys1[i] = verts(i, j)->y();
            }
            plt.plot(xs1, ys1, {.linewidth=1.0, .color="black"});
        }

        std::vector<double> xs2(block->size2() + 1);
        std::vector<double> ys2(block->size2() + 1);
        for (int i = 0; i <= block->size1(); ++i) {
            for (int j = 0; j <= block->size2(); ++j) {
                xs2[j] = verts(i, j)->x();
                ys2[j] = verts(i, j)->y();
            }
            plt.plot(xs2, ys2, {.linewidth=1.0, .color="black"});
        }

        for (int i = 0; i <= block->size1(); ++i) {
            for (int j = 0; j <= block->size2(); ++j) {
                // plt.text(verts(i, j)->x(), verts(i, j)->y(), std::to_string(verts(i,j)->degree()));
            }
        }

        // Связи с соседями
        for (int i = 0; i <= block->size1(); ++i) {
            for (int j = 0; j <= block->size2(); ++j) {
                Vector3d v1 = verts(i, j)->v1;
                for (auto edge: verts(i, j)->adjacent()) {
                    Vector3d v2 = edge.neib->v1;
                    Vector3d dr = 0.25 * (v2 - v1);

                    //plt.arrow(v1.x(), v1.y(), dr.x(), dr.y());

                    Vector3d vc = 0.5 * (v1 + v2);
                    auto ptr1 = reinterpret_cast<uintptr_t>(edge.ratio1);
                    auto ptr2 = reinterpret_cast<uintptr_t>(edge.ratio2);
                    if (ptr1 > ptr2) std::swap(ptr1, ptr2);
                    //plt.text(vc.x(), vc.y(), std::format("({}, {})", ptr1 % 97, ptr2 % 97), {.ha = "center", .va = "center"});
                }
            }
        }
    }

    plt.tight_layout();
    plt.show();
}

void BlockStructured::set_accuracy(double eps) {
    m_epsilon = std::max(0.0, std::min(eps, 1.0));
}

void BlockStructured::set_verbose(bool v) {
    m_verbose = v;
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

void BlockStructured::check_consistency() const {
    for (auto& block: m_blocks) {
        if (block->size1() < 1 || block->size2() < 1) {
            throw std::runtime_error("Zero size of some block");
        }

        for (Side2D side: Side2D::items()) {
            if (block->boundary(side)) {
                continue;
            }

            auto neib = block->adjacent_block(side);
            if (neib) {
                Side2D twin = block->neib_face(side);

                int size1 = block->size(side);
                int size2 = neib->size(twin);

                if (size1 != size2) {
                    std::cout << "side : " << side << "\n";
                    throw std::runtime_error("Sizes mismatch");
                }
            }
            else {
                throw std::runtime_error("Boundary or neighbor should be defined");
            }
        }
    }
}

void BlockStructured::initialize() const {
    check_consistency();

    for (const auto &block: m_blocks) {
        block->create_vertices();
    }

    for (auto& block: m_blocks) {
        block->link_vertices();
    }

    int counter = 0;
    double error = 1.0;
    if (m_verbose) {
        std::cout << "smoothing:\n";
    }

    // Фиксируем часть вершин по критерию
    for (auto& block: m_blocks) {
        for (int i = 0; i < block->m_vertices.size1(); ++i) {
            for (int j = 0; j < block->m_vertices.size2(); ++j) {
                auto node = block->m_vertices(i, j);
                if (m_fixed(node->v1)) { node->fix(); }
            }
        }
    }

    while (error > m_epsilon) {
        if (m_verbose && counter % 100 == 0) {
            std::cout << std::scientific << std::setprecision(2);
            std::cout << "\tstep " << std::setw(8) << counter << "\t\teps: " << error << "\n";
        }
        error = 0.0;
        for (auto &block: m_blocks) {
            double err = block->smooth();
            error = std::max(error, err);
        }
        for (auto &block: m_blocks) {
            block->update();
        }
        ++counter;
    }
}

void BlockStructured::optimize() {
    int N1 = 10;
    int N2 = 10;
    for (const auto& block: m_blocks) {
        block->set_size(0, N1);
        block->set_size(1, N2);
    }

    check_consistency();

    for (const auto &block: m_blocks) {
        block->create_vertices();
    }

    for (const auto& block: m_blocks) {
        block->link_vertices();
    }
    plot();

    //m_blocks[0]->m_modulus = 1.0;
    //m_blocks[1]->m_modulus = 1.0;

    int counter = 0;
    double error = 1.0;

    for (int i = 0; i < 5000; ++i) {
        if (m_verbose && counter % 100 == 0) {
            std::cout << std::scientific << std::setprecision(2);
            std::cout << "\tstep " << std::setw(8) << counter << "\t\teps: " << error << "\n";
        }
        error = 0.0;
        for (const auto &block: m_blocks) {
            double err = block->smooth();
            error = std::max(error, err);
        }
        for (const auto &block: m_blocks) {
            block->update();
            block->update_modulus();
        }
        ++counter;
    }
    plot();
}

void BlockStructured::optimize2() {
    // Очистим блоки
    for (auto& block: m_blocks) {
        block->set_size(0, 1);
        block->set_size(1, 1);
        block->m_vertices = {};
    }

    double max_modulus = 1.0;
    Block::Ptr bad_block;
    int axis = 0;
    for (auto& block: m_blocks) {
        if (block->modulus() > max_modulus) {
            max_modulus = block->modulus();
            bad_block = block;
            axis = 1;
        }
        if (1.0 / block->modulus() > max_modulus) {
            max_modulus = 1.0 / block->modulus();
            bad_block = block;
            axis = 0;
        }
    }
    std::cout << "max_modulus: " << max_modulus << "; bad block: " << bad_block->index << "\n";

    bad_block->set_size(axis, 10);

    for (int i = 0; i < 3; ++i) {
        for (auto block: m_blocks) {
            std::cout << "block sizes: " << block->size1() << " " << block->size2() << "\n";
            if (block->size1() == 1 && block->size2() > 1) {
                block->set_size(0, std::max(static_cast<int>(std::round(block->size2() * block->modulus())), 1));
                std::cout << "HERE " << block->size2() << " " <<  block->modulus() << " " << block->size1() << "\n";
            }
            if (block->size1() > 1 && block->size2() == 1) {
                block->set_size(1, std::max(static_cast<int>(std::round(block->size1() / block->modulus())), 1));
                std::cout << "HERE " << block->size1() << " " <<  block->modulus() << " " << block->size2() << "\n";
            }
        }
    }

    check_consistency();

    for (const auto &block: m_blocks) {
        block->create_vertices();
    }

    for (const auto& block: m_blocks) {
        block->link_vertices();
    }
    plot();
    //return;

    //m_blocks[0]->m_modulus = 1.0;
    //m_blocks[1]->m_modulus = 1.0;

    int counter = 0;
    double error = 1.0;

    for (int i = 0; i < 5000; ++i) {
        if (m_verbose && counter % 100 == 0) {
            std::cout << std::scientific << std::setprecision(2);
            std::cout << "\tstep " << std::setw(8) << counter << "\t\teps: " << error << "\n";
        }
        error = 0.0;
        for (const auto &block: m_blocks) {
            double err = block->smooth();
            error = std::max(error, err);
        }
        for (const auto &block: m_blocks) {
            block->update();
            block->update_modulus();
        }
        ++counter;
    }
    plot();
}

void BlockStructured::optimize3() {
    // Очистим блоки
    for (auto& block: m_blocks) {
        block->set_size(0, 1);
        block->set_size(1, 1);
        block->m_vertices = {};
    }

    double max_modulus = 1.0;
    Block::Ptr bad_block;
    int axis = 0;
    for (auto& block: m_blocks) {
        if (block->modulus() > max_modulus) {
            max_modulus = block->modulus();
            bad_block = block;
            axis = 1;
        }
        if (1.0 / block->modulus() > max_modulus) {
            max_modulus = 1.0 / block->modulus();
            bad_block = block;
            axis = 0;
        }
    }
    std::cout << "max_modulus: " << max_modulus << "; bad block: " << bad_block->index << "\n";

    bad_block->set_size(axis, 10);

    for (int i = 0; i < 3; ++i) {
        for (auto block: m_blocks) {
            std::cout << "block sizes: " << block->size1() << " " << block->size2() << "\n";
            if (block->size1() == 1 && block->size2() > 1) {
                block->set_size(0, std::max(static_cast<int>(std::round(block->size2() * block->modulus())), 1));
                std::cout << "HERE " << block->size2() << " " <<  block->modulus() << " " << block->size1() << "\n";
            }
            if (block->size1() > 1 && block->size2() == 1) {
                block->set_size(1, std::max(static_cast<int>(std::round(block->size1() / block->modulus())), 1));
                std::cout << "HERE " << block->size1() << " " <<  block->modulus() << " " << block->size2() << "\n";
            }
        }
    }

    check_consistency();

    for (const auto &block: m_blocks) {
        block->create_vertices();
    }

    for (const auto& block: m_blocks) {
        block->link_vertices();
    }
    plot();
    //return;

    //m_blocks[0]->m_modulus = 1.0;
    //m_blocks[1]->m_modulus = 1.0;

    int counter = 0;
    double error = 1.0;

    for (auto block: m_blocks) {
        block->m_ratio = {1.0, 1.0};
    }

    for (int i = 0; i < 5000; ++i) {
        if (m_verbose && counter % 100 == 0) {
            std::cout << std::scientific << std::setprecision(2);
            std::cout << "\tstep " << std::setw(8) << counter << "\t\teps: " << error << "\n";
        }
        error = 0.0;
        for (const auto &block: m_blocks) {
            double err = block->smooth();
            error = std::max(error, err);
        }
        for (const auto &block: m_blocks) {
            block->update();
            // block->update_modulus();
        }
        ++counter;
    }
    plot();
}

Grid BlockStructured::make() {
    initialize();

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

} // namespace zephyr::geom::generator
