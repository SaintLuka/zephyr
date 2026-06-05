#include <format>
#include <zephyr/geom/box.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/block_structured.h>
#include <zephyr/utils/pyplot.h>

namespace zephyr::geom::generator {

inline std::tuple<double, double> figsize(Box bbox) {
    constexpr double max_width{10.0};
    constexpr double max_height{6.0};

    double aspect = bbox.sizes()[0] / bbox.sizes()[1];

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

void BlockStructured::plot_layout() const {
    // Собираем уникальные вершины
    Box bbox = Box::Empty(2);
    std::set<BaseNode::Ptr> unique_nodes;
    for (const auto& block: m_blocks) {
        for (const auto& node: block->base_nodes()) {
            unique_nodes.insert(node);
            bbox.capture(node->pos());
        }
    }

    utils::pyplot plt;
    plt.figure({.figsize = figsize(bbox), .dpi = 175});
    plt.set_aspect_equal();

    // Очертания блоков
    for (const auto& block: m_blocks) {
        std::vector<double> bxs, bys;
        for (int i: {0, 1, 3, 2, 0}) {
            bxs.push_back(block->base_node(i)->x());
            bys.push_back(block->base_node(i)->y());
        }
        plt.plot(bxs, bys, {.linestyle="solid", .color="black", .marker="."});
    }

    // Построить узлы
    std::vector<double> nodes_x, nodes_y;
    for (const auto& node: unique_nodes) {
        nodes_x.push_back(node->x());
        nodes_y.push_back(node->y());
    }
    plt.plot(nodes_x, nodes_y, {.linestyle="none", .marker="o"});

    // Центры блоков
    for (int ib = 0; ib < m_blocks.size(); ++ib) {
        Vector3d bc = m_blocks[ib]->center();
        plt.text(bc.x(), bc.y(), std::format("$B_{{{}}}$", ib), {.ha = "center", .va = "center"});
    }

    plt.tight_layout();
    plt.show();
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

    // Центры блоков
    /*
    for (int ib = 0; ib < m_blocks.size(); ++ib) {
        Vector3d bc = m_blocks[ib]->center();
        plt.text(bc.x(), bc.y(), std::format("$B_{{{}}}$", ib), {.ha = "center", .va = "center"});
    }
    */

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

    // Центры блоков
    /*
    for (int ib = 0; ib < all_vertices.size(); ++ib) {
        int i = all_vertices[ib].size(Axis::X) / 2;
        int j = all_vertices[ib].size(Axis::Y) / 2;
        Vector3d bc = all_vertices[ib](i, j)->pos;
        plt.text(bc.x(), bc.y(), std::format("$B_{{{}}}$", ib), {.ha = "center", .va = "center"});
    }
    */

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

} // namespace zephyr::geom::generator
