#include <array>
#include <iomanip>

#include <zephyr/mesh/storage.h>

#include <zephyr/geom/cell.h>

#include <zephyr/mesh/generator/block.h>
#include <zephyr/mesh/generator/block_structured.h>
#include <zephyr/mesh/generator/vertex.h>
#include <zephyr/mesh/generator/curve/curve.h>


namespace zephyr { namespace mesh { namespace generator {

BlockStructured::BlockStructured(int n_blocks)
    : Generator("block_structured") {

    m_blocks.reserve(n_blocks);
    for (int index = 0; index < n_blocks; ++index) {
        m_blocks.push_back(Block(index));
    }

    m_epsilon = 1.0e-3;
    m_verbose = false;
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
    int count = 0;
    for (auto& block: m_blocks) {
        count += block.size1() * block.size2();
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

VerticesList BlockStructured::get_vertices() {
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

    VerticesList vlist;
    for (auto &block: m_blocks) {
        block.add_vertices(vlist);
    }
    return vlist;
}

BoundaryFlags BlockStructured::get_boundary_flags() {
    BoundaryFlags flist;
    for (auto &block: m_blocks) {
        block.add_boundary_flags(flist);
    }
    return flist;
}

void BlockStructured::initialize(Storage& cells) {
    using zephyr::geom::Cell;
    using zephyr::geom::ShortList2D;

    auto vlist = get_vertices();
    auto flist = get_boundary_flags();

    cells.resize(vlist.size());
    for (int n = 0; n < cells.size(); ++n) {
        ShortList2D verts = vlist[n];

        Cell g_cell(verts);

        for (int s = 0; s < 4; ++s) {
            g_cell.faces[s].boundary = flist[n][s];
        }

        cells[n].geom() = g_cell;
    }
}

} // namespace generator
} // namespace mesh
} // namespace zephyr
