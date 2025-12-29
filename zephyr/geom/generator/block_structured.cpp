#include <array>
#include <algorithm>
#include <iomanip>

#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>
#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/block_structured.h>
#include <zephyr/utils/pyplot.h>


namespace zephyr::geom::generator {

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

void BlockStructured::link() {
    if (m_stage != EDITABLE) {
        throw std::runtime_error("BlockStructured::link() is available only at edit stage");
    }

    // Удаляем нулевые блоки
    const auto it = std::remove(m_blocks.begin(), m_blocks.end(), nullptr);
    m_blocks.erase(it, m_blocks.end());

    // Собираем уникальные вершины, чистим ссылки на блоки
    std::set<BaseVertex::Ptr> nodes;
    for (const auto& block: m_blocks) {
        for (int i = 0; i < 4; ++i) {
            auto v = block->base_vertex(i);
            v->clear_adjacent_blocks();
            nodes.insert(v);
        }
    }

    // Для всех вершин выставляем смежные блоки
    for (auto& block: m_blocks) {
        for (int i = 0; i < 4; ++i) {
            auto v = block->base_vertex(i);
            v->add_adjacent_block(block);
        }
    }

    // Для каждой вершины проходим по смежным блокам
    for (auto &bv: nodes) {
        auto &adj = bv->adjacent_blocks();
        for (auto it_b1 = adj.begin(); it_b1 != adj.end(); ++it_b1) {
            if (it_b1->expired()) continue;

            auto it_b2 = it_b1; ++it_b2;
            for (; it_b2 != adj.end(); ++it_b2) {
                if (it_b2->expired()) continue;

                it_b1->lock()->link(it_b2->lock());
            }
        }
    }

    m_stage = LINKED;
}

void BlockStructured::plot() const {
    // Собираем уникальные вершины
    std::set<BaseVertex::Ptr> nodes;
    for (const auto& block: m_blocks) {
        for (int i = 0; i < 4; ++i) {
            nodes.insert(block->base_vertex(i));
        }
    }

    std::vector<double> nxs; nxs.reserve(nodes.size());
    std::vector<double> nys; nys.reserve(nodes.size());
    for (auto node: nodes) {
        nxs.push_back(node->v().x());
        nys.push_back(node->v().y());
    }

    utils::pyplot plt;
    plt.figure({.dpi=170});

    plt.set_aspect_equal();
    //plt.plot(nxs, nys, {.linestyle="none", .marker="o"});

    std::vector<double> bxs(5);
    std::vector<double> bys(5);
    for (int ib = 0; ib < m_blocks.size(); ++ib) {
        auto& block = m_blocks[ib];

        /*
        for (int iv = 0; iv <= 4; ++iv) {
            bxs[iv] = block->base_vertex(iv % 4)->v().x();
            bys[iv] = block->base_vertex(iv % 4)->v().y();
        }
        double xc = 0.25 * std::accumulate(bxs.begin(), bxs.end() - 1, 0.0);
        double yc = 0.25 * std::accumulate(bys.begin(), bys.end() - 1, 0.0);

        plt.plot(bxs, bys, {.linestyle="solid", .marker="."});
        plt.text(xc, yc, "B" + std::to_string(ib));
        */

        for (int iv = 0; iv < 4; ++iv) {
            double x = block->base_vertex(iv % 4)->v().x();
            double y = block->base_vertex(iv % 4)->v().y();
            //plt.text(x, y, std::to_string(block->base_vertex(iv)->adjacent_blocks().size()));
        }

        const auto& verts = block->m_vertices;

        std::vector<double> xs1(block->size1() + 1);
        std::vector<double> ys1(block->size1() + 1);
        for (int j = 0; j <= block->size2(); ++j) {
            for (int i = 0; i <= block->size1(); ++i) {
                xs1[i] = verts[i][j]->v1.x();
                ys1[i] = verts[i][j]->v1.y();
            }
            plt.plot(xs1, ys1, {.linewidth=1.0, .color="black"});
        }

        std::vector<double> xs2(block->size2() + 1);
        std::vector<double> ys2(block->size2() + 1);
        for (int i = 0; i <= block->size1(); ++i) {
            for (int j = 0; j <= block->size2(); ++j) {
                xs2[j] = verts[i][j]->v1.x();
                ys2[j] = verts[i][j]->v1.y();
            }
            plt.plot(xs2, ys2, {.linewidth=1.0, .color="black"});
        }

        for (int i = 0; i <= block->size1(); ++i) {
            for (int j = 0; j <= block->size2(); ++j) {
                //plt.text(verts[i][j]->v1.x(), verts[i][j]->v2.y(), std::to_string(verts[i][j]->n_adjacent()));
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

        for (int f1 = 0; f1 < 4; ++f1) {
            if (block->boundary(f1)) {
                continue;
            }

            Block::Ptr neib = block->adjacent_block(f1);
            std::cout << "side : " << f1 << "\n";

            if (neib) {
                int f2 = block->neib_face(f1);

                int size1 = block->size(f1);
                int size2 = neib->size(f2);

                if (size1 != size2) {
                    throw std::runtime_error("Sizes mismatch");
                }
            }
            else {
                throw std::runtime_error("Boundary or neighbor should be defined");
            }
        }
    }
}

void BlockStructured::initialize() {
    check_consistency();

    for (auto &block: m_blocks) {
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
        for (auto &row: block->m_vertices) {
            for (auto &node: row) {
                if (m_fixed(node->v1)) {
                    node->fix();
                }
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
    int N1 = 20;
    int N2 = 20;
    for (auto& block: m_blocks) {
        block->set_size_to_face(0, N1);
        block->set_size_to_face(1, N2);
    }

    check_consistency();

    for (auto &block: m_blocks) {
        block->create_vertices();
    }

    for (auto& block: m_blocks) {
        block->link_vertices();
    }

    m_blocks[0]->m_modulus = N2 / double(N1);
    //plot();

    int counter = 0;
    double error = 1.0;

    for (int i = 0; i < 5500; ++i) {
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

            // Update modulus
            double M = block->calc_modulus();
            block->m_modulus = block->size2() / (M * block->size1());
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
