#include <array>
#include <iomanip>

#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>
#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/block_structured.h>


namespace zephyr::geom::generator {

BlockStructured::BlockStructured(int n_blocks)
    : Generator("block_structured") {

    m_blocks.reserve(n_blocks);
    for (int index = 0; index < n_blocks; ++index) {
        m_blocks.push_back(Block(index));
    }

    m_epsilon = 1.0e-3;
    m_verbose = false;

    m_fixed = [](const Vector3d&) -> bool {
        return false;
    };
}

Block &BlockStructured::operator[](int idx) {
    if (idx >= m_blocks.size()) {
        throw std::out_of_range("BlockStructured::operator[]");
    }
    return m_blocks[idx];
}

const Block &BlockStructured::operator[](int idx) const {
    if (idx >= m_blocks.size()) {
        throw std::out_of_range("BlockStructured::operator[]");
    }
    return m_blocks[idx];
}

void BlockStructured::link() {
    std::set<BaseVertex::Ptr> verts;

    for (auto &block: m_blocks) {
        for (int i = 0; i < 4; ++i) {
            BaseVertex::Ptr v = block.base_vertex(i);
            v->add_adjacent_block(&block);
            verts.insert(v);
        }
    }

    for (auto &bv: verts) {
        auto &adj = bv->adjacent_blocks();
        for (int i = 0; i < adj.size(); ++i) {
            for (int j = i + 1; j < adj.size(); ++j) {
                adj[i]->link(adj[j]);
            }
        }
    }
}

void BlockStructured::set_accuracy(double eps) {
    m_epsilon = std::max(0.0, std::min(eps, 1.0));
}

void BlockStructured::set_verbose(bool v) {
    m_verbose = v;
}

int BlockStructured::size() const {
    return calc_cells();
}

int BlockStructured::calc_cells() const {
    int count = 0;
    for (auto& block: m_blocks) {
        count += block.size1() * block.size2();
    }
    return count;
}

int BlockStructured::calc_nodes() const {
    int count = 0;
    for (auto& block: m_blocks) {
        count += (block.size1() + 1) * (block.size2() + 1);
    }
    return count;
}

void BlockStructured::check_consistency() const {
    for (auto& block: m_blocks) {
        if (block.size1() < 1 || block.size2() < 1) {
            throw std::runtime_error("Zero size of some block");
        }

        for (int f1 = 0; f1 < 4; ++f1) {
            if (block.boundary(f1)) {
                continue;
            }

            Block* neib = block.adjacent_block(f1);

            if (neib) {
                int f2 = block.neib_face(f1);

                int size1 = block.size(f1);
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
        block.create_vertices();
    }

    for (auto& block: m_blocks) {
        block.link_vertices();
    }

    int counter = 0;
    double error = 1.0;
    if (m_verbose) {
        std::cout << "smoothing:\n";
    }

    // Фиксируем часть вершин по критерию
    for (auto& block: m_blocks) {
        for (auto &row: block.m_vertices) {
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
            double err = block.smooth();
            error = std::max(error, err);
        }
        for (auto &block: m_blocks) {
            block.update();
        }
        ++counter;
    }
}

Grid BlockStructured::make() {
    initialize();

    Grid grid;
    grid.reserve_cells(calc_cells());
    grid.reserve_nodes(calc_nodes());

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

    return grid;
}

} // namespace zephyr::geom::generator
