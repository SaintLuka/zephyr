#include <algorithm>
#include <array>
#include <map>
#include <set>

#include <zephyr/geom/maps.h>

#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>


namespace zephyr::geom::generator {

inline int next(int idx) {
    return (idx + 1) % 4;
}

inline int prev(int idx) {
    return (idx + 3) % 4;
}

inline int opp(int idx) {
    return (idx + 2) % 4;
}

Block::Block(int index)
    : m_index(index) {

    m_size1 = 0;
    m_size2 = 0;
    m_boundaries = {nullptr, nullptr, nullptr, nullptr};
    m_base_vertices = {nullptr, nullptr, nullptr, nullptr};
    m_adjacent_blocks = {nullptr, nullptr, nullptr, nullptr};
    m_rotations = {-1, -1, -1, -1};
}

int Block::index() const {
    return m_index;
}

int Block::size(BaseVertex::Ref bv1, BaseVertex::Ref bv2) const {
    return face_index(bv1, bv2) % 2 ? m_size2 : m_size1;
}

int Block::size1() const {
    return m_size1;
}

int Block::size2() const {
    return m_size2;
}

Block &Block::operator=(std::initializer_list<BaseVertex::Ptr> vertices) {
    if (vertices.size() != 4) {
        throw std::runtime_error("Wrong number of base vertices");
    }

    std::vector<std::pair<double, BaseVertex::Ptr>> vec_verts;
    for (auto &v: vertices) {
        vec_verts.emplace_back(std::make_pair(0.0, v));
    }

    Vector3d v0 = {0.0, 0.0, 0.0};
    for (int i = 0; i < 4; ++i) {
        v0 += vec_verts[i].second->v();
    }
    v0 /= 4;

    for (int i = 0; i < 4; ++i) {
        vec_verts[i].first =
            std::atan2(vec_verts[i].second->v().y() - v0.y(),
                       vec_verts[i].second->v().x() - v0.x());

        if (vec_verts[i].first < vec_verts[0].first) {
            vec_verts[i].first += 2.0 * M_PI;
        }
    }

    std::sort(vec_verts.begin(), vec_verts.end(),
              [](const std::pair<double, BaseVertex::Ptr> &p1, const
              std::pair<double, BaseVertex::Ptr> &p2) {
                  return p1.first < p2.first;
              });

    for (int i = 0; i < 4; ++i) {
        m_base_vertices[i] = vec_verts[i].second;
    }

    return *this;
}

bool Block::is_boundary(BaseVertex::Ref v) const {
    for (auto block: v->adjacent_blocks()) {
        int v_idx = block->vertex_index(v);
        if (block->boundary(v_idx) != nullptr ||
            block->boundary(prev(v_idx)) != nullptr) {
            return true;
        }
    }
    return false;
}

Curve::Ref Block::boundary(int f_idx) const {
    return m_boundaries[f_idx];
}

Block* Block::adjacent_block(int f_idx) const {
    return m_adjacent_blocks[f_idx];
}

BaseVertex::Ptr Block::base_vertex(int v_idx) const {
    return m_base_vertices[v_idx];
}

int Block::vertex_index(BaseVertex::Ref v) const {
    int idx = 0;
    while (idx < 4 && m_base_vertices[idx] != v) {
        ++idx;
    }
    if (idx >= 4) {
        throw std::runtime_error("Can't find base vertex");
    }
    return idx;
}

int Block::face_index(BaseVertex::Ref v1, BaseVertex::Ref v2) const {
    int idx = vertex_index(v1);

    if (v2 == m_base_vertices[next(idx)]) {
        return idx;
    }
    else if (v2 == m_base_vertices[prev(idx)]) {
        return prev(idx);
    }
    else {
        throw std::runtime_error("Can't find face between two vertices");
    }
}

void Block::set_boundary(BaseVertex::Ref v1, BaseVertex::Ref v2, Curve::Ref curve) {
    m_boundaries[face_index(v1, v2)] = curve;
}

BsVertex::Ptr Block::vertex(int i, int j) const {
    assert(i <= m_size1);
    assert(j <= m_size2);

    return m_vertices[i][j];
}

int Block::size(int f_idx) const {
    return f_idx % 2 ? m_size2 : m_size1;
}

void Block::set_size(BaseVertex::Ref v1, BaseVertex::Ref v2, int N) {
    int k = face_index(v1, v2);
    if (k % 2) {
        m_size2 = N;
    } else {
        m_size1 = N;
    }

    for (int l: {k, opp(k)}) {
        if (m_adjacent_blocks[l]) {
            BaseVertex::Ref BV1 = m_base_vertices[l];
            BaseVertex::Ref BV2 = m_base_vertices[next(l)];
            if (m_adjacent_blocks[l]->size(BV1, BV2) != N) {
                m_adjacent_blocks[l]->set_size(BV1, BV2, N);
            }
        }
    }
}

BsVertex::Ptr &Block::corner_vertex(int v_idx) {
    switch (v_idx) {
        case 0:
            return m_vertices[0][0];
        case 1:
            return m_vertices[m_size1][0];
        case 2:
            return m_vertices[m_size1][m_size2];
        default:
            return m_vertices[0][m_size2];
    }
}

BsVertex::Ptr &Block::boundary_vertex(int f_idx, int idx) {
    switch (f_idx) {
        case 0:
            return m_vertices[idx][0];
        case 1:
            return m_vertices[m_size1][idx];
        case 2:
            return m_vertices[m_size1 - idx][m_size2];
        default:
            return m_vertices[0][m_size2 - idx];
    }
}

BsVertex::Ptr &Block::preboundary_vertex(int f_idx, int idx) {
    switch (f_idx) {
        case 0:
            return m_vertices[idx][1];
        case 1:
            return m_vertices[m_size1 - 1][idx];
        case 2:
            return m_vertices[m_size1 - idx][m_size2 - 1];
        default:
            return m_vertices[1][m_size2 - idx];
    }
};

int Block::neib_face(int f_idx) const {
    return (f_idx - m_rotations[f_idx] + 4) % 4;
}

void Block::link(Block* block) {
    for (int f1 = 0; f1 < 4; ++f1) {
        BaseVertex::Ptr a1 = base_vertex(f1);
        BaseVertex::Ptr b1 = base_vertex(next(f1));

        for (int f2 = 0; f2 < 4; ++f2) {
            BaseVertex::Ptr a2 = block->base_vertex(f2);
            BaseVertex::Ptr b2 = block->base_vertex(next(f2));

            if (a1 == b2 && b1 == a2) {
                m_adjacent_blocks[f1] = block;
                block->m_adjacent_blocks[f2] = this;

                m_rotations[f1] = (f1 - f2 + 4) % 4;
                block->m_rotations[f2] = (f2 - f1 + 4) % 4;
            }
        }
    }
}

void Block::create_vertices() {
    m_vertices.resize(m_size1 + 1, std::vector<BsVertex::Ptr>(m_size2 + 1, nullptr));

    using zephyr::geom::SqQuad;

    Vector3d v0 = base_vertex(0)->v();
    Vector3d v1 = base_vertex(1)->v();
    Vector3d v2 = base_vertex(3)->v();
    Vector3d v3 = base_vertex(2)->v();
    Vector3d vL = (v0 + v2) / 2.0;
    Vector3d vR = (v1 + v3) / 2.0;
    Vector3d vB = (v0 + v1) / 2.0;
    Vector3d vT = (v2 + v3) / 2.0;

    if (m_boundaries[3]) {
        vL = m_boundaries[3]->projection(vL);
    }
    if (m_boundaries[1]) {
        vR = m_boundaries[1]->projection(vR);
    }
    if (m_boundaries[0]) {
        vB = m_boundaries[0]->projection(vB);
    }
    if (m_boundaries[2]) {
        vT = m_boundaries[2]->projection(vT);
    }

    Vector3d vC = (vL + vR + vB + vT) / 4.0;
    SqQuad quad = {
            v0, vB, v1,
            vL, vC, vR,
            v2, vT, v3
    };

    for (int i = 0; i <= m_size1; ++i) {
        double x = (2.0 * i - m_size1) / m_size1;

        m_vertices[i][0]       = BsVertex::create(quad(x, -1.0));
        m_vertices[i][m_size2] = BsVertex::create(quad(x, +1.0));
    }

    for (int j = 0; j <= m_size2; ++j) {
        double y = (2.0 * j - m_size2) / m_size2;

        m_vertices[0][j]       = BsVertex::create(quad(-1.0, y));
        m_vertices[m_size1][j] = BsVertex::create(quad(+1.0, y));
    }

    for (int k = 0; k < 0; ++k) {
        for (int f_idx = 0; f_idx < 4; ++f_idx) {
            if (!boundary(f_idx)) {
                continue;
            }

            for (int idx = 1; idx < size(f_idx); ++idx) {
                Vector3d v = (boundary_vertex(f_idx, idx - 1)->v1 +
                    boundary_vertex(f_idx, idx + 1)->v1) / 2.0;
                v = boundary(f_idx)->projection(v);

                boundary_vertex(f_idx, idx)->v2 = v;
            }
        }

        for (int f_idx = 0; f_idx < 4; ++f_idx) {
            if (!boundary(f_idx)) {
                continue;
            }

            for (int idx = 1; idx < size(f_idx); ++idx) {
                boundary_vertex(f_idx, idx)->v1 = boundary_vertex(f_idx, idx)->v2;
            }
        }
    }

    for (int i = 1; i < m_size1; ++i) {
        for (int j = 1; j < m_size2; ++j) {
            double x = (2.0 * i - m_size1) / m_size1;
            double y = (2.0 * j - m_size2) / m_size2;

            m_vertices[i][j] = BsVertex::create(quad(x, y));
        }
    }

    for (int i = 0; i < 4; ++i) {
        auto neib = m_adjacent_blocks[i];

        if (!neib) {
            continue;
        }

        // Сосед есть и уже с вершинами
        if (neib->m_vertices.empty()) {
            continue;
        }

        int j = neib_face(i);

        int N = size(i);
        for (int k = 0; k <= N; ++k) {
            boundary_vertex(i, k) = neib->boundary_vertex(j, N - k);
        }
    }
}

void Block::link_vertices() {
    // Базисные вершины
    for (int v_idx = 0; v_idx < 4; ++v_idx) {
        auto bv = m_base_vertices[v_idx];

        /// Фиксированная точка
        if (bv->is_fixed()) {
            continue;
        }

        // Угловая вершина области, должна быть неподвижной
        if (bv->degree() < 2) {
            continue;
        }

        // Угловая вершина области по другому критерию
        if (m_boundaries[v_idx] && m_boundaries[prev(v_idx)]) {
            continue;
        }

        Curve *boundary = nullptr;
        std::vector<BsVertex::Ptr> adj_vertices;

        // Вершина на границе
        if (is_boundary(bv)) {
            int f_idx = v_idx;
            Block* block = this;

            while (block) {
                adj_vertices.insert(
                    adj_vertices.begin(),
                    block->boundary_vertex(f_idx, 1)
                );

                boundary = block->boundary(f_idx).get();

                Block *next_block = block->adjacent_block(f_idx);
                f_idx = next(block->neib_face(f_idx));
                block = next_block;

                if (adj_vertices.size() > bv->degree() + 1) {
                    throw std::runtime_error("Infinite loop #1");
                }
            }

            f_idx = prev(v_idx);
            block = this;

            while (block) {
                adj_vertices.push_back(
                    block->boundary_vertex(f_idx, block->size(f_idx) - 1)
                );

                boundary = block->boundary(f_idx).get();

                Block *next_block = block->adjacent_block(f_idx);
                f_idx = prev(block->neib_face(f_idx));
                block = next_block;

                if (adj_vertices.size() > bv->degree() + 1) {
                    throw std::runtime_error("Infinite loop #2");
                }
            }

            if (adj_vertices.size() != bv->degree() + 1) {
                throw std::runtime_error("Boundary error");
            }

        } else {
            for (auto adj: bv->adjacent_blocks()) {
                auto v_idx2 = adj->vertex_index(bv);
                adj_vertices.emplace_back(adj->boundary_vertex(v_idx2, 1));
            }
        }

        corner_vertex(v_idx)->set_boundary(boundary);
        corner_vertex(v_idx)->set_adjacent_vertices(adj_vertices);
    }

    // Вершины на границе
    for (int f_idx = 0; f_idx < 4; ++f_idx) {
        if (!boundary(f_idx) && !adjacent_block(f_idx)) {
            throw std::runtime_error("AmrFace is not boundary and has no neigbor");
        }

        // Число ячеек по грани
        int N = size(f_idx);

        if (boundary(f_idx)) {

            for (int idx = 1; idx < N; ++idx) {
                boundary_vertex(f_idx, idx)->set_adjacent_vertices(
                    {
                        boundary_vertex(f_idx, idx + 1),
                        preboundary_vertex(f_idx, idx),
                        boundary_vertex(f_idx, idx - 1)
                    }
                );

                boundary_vertex(f_idx, idx)->set_boundary(
                    m_boundaries[f_idx].get()
                );
            }

        } else {
            auto neib = m_adjacent_blocks[f_idx];
            int fn_idx = neib_face(f_idx);

            for (int idx = 1; idx < N; ++idx) {
                auto K0 = boundary_vertex(f_idx, idx);
                auto K1 = boundary_vertex(f_idx, idx - 1);
                auto K2 = boundary_vertex(f_idx, idx + 1);
                auto K3 = preboundary_vertex(f_idx, idx);
                auto K4 = neib->preboundary_vertex(fn_idx, N - idx);

                boundary_vertex(f_idx, idx)->set_adjacent_vertices(
                    {
                        boundary_vertex(f_idx, idx - 1),
                        boundary_vertex(f_idx, idx + 1),
                        preboundary_vertex(f_idx, idx),
                        neib->preboundary_vertex(fn_idx, N - idx)
                    });
            }
        }
    }

    // Внутренние вершины
    for (int i = 1; i < m_size1; ++i) {
        for (int j = 1; j < m_size2; ++j) {
            m_vertices[i][j]->set_adjacent_vertices(
                {
                    m_vertices[i + 1][j],
                    m_vertices[i - 1][j],
                    m_vertices[i][j - 1],
                    m_vertices[i][j + 1]
                }
            );
        }
    }
}

double Block::smooth() {
    double err = 0.0;

    for (auto& row: m_vertices) {
        for (auto& vertex: row) {
            Curve *boundary = vertex->boundary();
            auto &adjacent = vertex->adjacent_vertices();

            Vector3d avg = {0.0, 0.0, 0.0};

            if (adjacent.empty()) {
                // Фиксированная вершина
                avg = vertex->v1;
            } else if (boundary) {
                // Вершина на границе
                avg += adjacent[0]->v1 / 2.0;
                for (int n = 1; n < int(adjacent.size()) - 1; ++n) {
                    avg += boundary->projection(adjacent[n]->v1);
                }
                avg += adjacent.back()->v1 / 2.0;
                avg /= adjacent.size() - 1.0;

                avg = vertex->boundary()->projection(avg);
            } else {
                // Внутренняя вершина
                for (auto neib: adjacent) {
                    avg += neib->v1;
                }
                avg /= vertex->n_adjacent();
            }

            vertex->v2 = avg;

            double L = 0.0;
            for (auto adj: adjacent) {
                L = std::max(L, (vertex->v1 - adj->v1).norm());
            }

            err = std::max(err, (vertex->v1 - vertex->v2).norm() / L);
        }
    }

    return err;
}

void Block::update() {
    for (auto &rows: m_vertices) {
        for (auto &vertex: rows) {
            vertex->v1 = vertex->v2;
        }
    }
}

} // namespace zephyr::geom::generator
