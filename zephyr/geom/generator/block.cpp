#include <algorithm>
#include <array>
#include <map>
#include <set>

#include <zephyr/geom/primitives/cube.h>

#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>

#include <zephyr/geom/intersection.h>
#include <zephyr/geom/primitives/triangle.h>

namespace zephyr::geom::generator {

constexpr int next(int idx) { return (idx + 1) % 4; }

constexpr int prev(int idx) { return (idx + 3) % 4; }

// Стороны вдоль первой оси (оси X)
constexpr bool axes1(Side2D side) { return side == Side2D::B || side == Side2D::T; }

constexpr Side2D next(Side2D side) {
    switch (side) {
        case Side2D::L: return Side2D::B;
        case Side2D::R: return Side2D::T;
        case Side2D::B: return Side2D::R;
        case Side2D::T: return Side2D::L;
        default: throw std::runtime_error("Invalid side");
    }
}

constexpr Side2D prev(Side2D side) {
    switch (side) {
        case Side2D::L: return Side2D::T;
        case Side2D::R: return Side2D::B;
        case Side2D::B: return Side2D::L;
        case Side2D::T: return Side2D::R;
        default: throw std::runtime_error("Invalid side");
    }
}

constexpr Side2D opp(Side2D side) {
    switch (side) {
        case Side2D::L: return Side2D::R;
        case Side2D::R: return Side2D::L;
        case Side2D::B: return Side2D::T;
        case Side2D::T: return Side2D::B;
        default: throw std::runtime_error("Invalid side");
    }
}

constexpr int node(Side2D side, int idx) {
    switch (side) {
        case Side2D::L: return std::array{0, 2}[idx];
        case Side2D::R: return std::array{1, 3}[idx];
        case Side2D::B: return std::array{0, 1}[idx];
        case Side2D::T: return std::array{2, 3}[idx];
        default: throw std::runtime_error("Invalid side");
    }
}

int Block::size(BaseNode::Ref bv1, BaseNode::Ref bv2) const {
    return axes1(get_side(bv1, bv2)) ? m_size1 : m_size2;
}

Block &Block::operator=(const std::array<BaseNode::Ptr, 4>& vertices) {
    if (vertices.size() != 4) {
        throw std::runtime_error("Wrong number of base vertices");
    }

    // Центр блока
    Vector3d vc = Vector3d::Zero();
    for (const auto& v: vertices) {
        vc += v->pos();
    }
    vc *= 0.25;

    // Пары (угол вершина)
    std::vector<std::pair<double, BaseNode::Ptr>> sorted(4);
    for (int i = 0; i < 4; ++i) {
        sorted[i].second = vertices[i];
        sorted[i].first = std::atan2(sorted[i].second->y() - vc.y(),
                                     sorted[i].second->x() - vc.x());
        if (sorted[i].first < sorted[0].first) {
            sorted[i].first += 2.0 * M_PI;
        }
    }

    // Сортировка по углу
    std::sort(sorted.begin() + 1, sorted.end(),
              [](const std::pair<double, BaseNode::Ptr> &p1,
                 const std::pair<double, BaseNode::Ptr> &p2) {
                  return p1.first < p2.first;
              });

    // Z-ordering
    m_base_nodes[0] = sorted[0].second;
    m_base_nodes[1] = sorted[1].second;
    m_base_nodes[2] = sorted[3].second;
    m_base_nodes[3] = sorted[2].second;
    return *this;
}

std::tuple<BaseNode::Ptr, BaseNode::Ptr> Block::adjacent_nodes(const BaseNode* v) const {
    if (!v) {
        throw std::runtime_error("Block::adjacent_nodes() error: Invalid BaseNode");
    }
    switch (node_index(v)) {
        case 0: return {base_node(1), base_node(2)};
        case 1: return {base_node(3), base_node(0)};
        case 2: return {base_node(0), base_node(3)};
        case 3: return {base_node(2), base_node(1)};
        default: throw std::runtime_error("Invalid index");
    }
}

constexpr std::tuple<Side2D, Side2D> node_adjacent_sides(int v_idx) {
    switch (v_idx) {
        case 0: return {Side2D::B, Side2D::L};
        case 1: return {Side2D::R, Side2D::B};
        case 2: return {Side2D::L, Side2D::T};
        case 3: return {Side2D::T, Side2D::R};
        default: throw std::runtime_error("Invalid index");
    }
}

std::tuple<Side2D, Side2D> Block::adjacent_sides(const BaseNode* v) const {
    if (!v) {
        throw std::runtime_error("Block::adjacent_sides() error: Invalid BaseNode");
    }
    return node_adjacent_sides(node_index(v));
}

std::tuple<Side2D, Side2D> Block::adjacent_sides(BaseNode::Ref v) const {
    return adjacent_sides(v.get());
}

Block::Ptr Block::adjacent_block(Side2D side) const {
    return m_adjacent_blocks[side].lock();
}

int Block::node_index(const BaseNode* v) const {
    for (int idx = 0; idx < 4; ++idx) {
        if (m_base_nodes[idx].get() == v) {
            return idx;
        }
    }
    throw std::runtime_error("Can't find base node");
}

bool Block::is_boundary(const BaseNode* v) const {
    auto[side1, side2] = adjacent_sides(v);
    return is_boundary(side1) || is_boundary(side2);
}

Side2D Block::get_side(BaseNode::Ref v1, BaseNode::Ref v2) const {
    // Глупый полный перебор, хотя тут всего по 4 узла/грани
    int idx_v1 = node_index(v1);
    int idx_v2 = node_index(v2);
    if (idx_v1 == idx_v2) {
        throw std::runtime_error("Invalid side, same nodes");
    }
    if (idx_v2 < idx_v1) {
        std::swap(idx_v1, idx_v2);
    }

    for (Side2D side: Side2D::items()) {
        int idx_w1 = node(side, 0);
        int idx_w2 = node(side, 1);
        if (idx_w2 < idx_w1) {
            std::swap(idx_w1, idx_w2);
        }

        if (idx_v1 == idx_w1 && idx_v2 == idx_w2) {
            return side;
        }
    }

    for (Side2D side: Side2D::items()) {
        int idx_w1 = node(side, 0);
        int idx_w2 = node(side, 1);
        if (idx_w2 < idx_w1) {
            std::swap(idx_w1, idx_w2);
        }

        if (idx_v1 == idx_w1 && idx_v2 == idx_w2) {
            return side;
        }
    }
    throw std::runtime_error("Can't find face between two vertices");
}

Vector3d Block::center() const {
    return 0.25 * (m_base_nodes[0]->pos() + m_base_nodes[1]->pos() + m_base_nodes[2]->pos() + m_base_nodes[3]->pos());
}

Vector3d Block::center(Side2D side) const {
    Vector3d v1 = base_node(node(side, 0))->pos();
    Vector3d v2 = base_node(node(side, 1))->pos();
    Vector3d vc = 0.5 * (v1 + v2);
    if (m_boundaries[side]) {
       vc = m_boundaries[side]->projection(vc);
    }
    return vc;
}

void Block::set_boundary(BaseNode::Ref v1, BaseNode::Ref v2, Curve::Ref curve) {
    auto side = get_side(v1, v2);
    m_boundaries[side] = curve;
    m_adjacent_blocks[side].reset();
}

BsVertex::Ptr Block::vertex(int i, int j) const {
    z_assert(i <= m_size1, "Block vertex: bad index");
    z_assert(j <= m_size2, "Block vertex: bad index");

    return m_vertices(i, j);
}

int Block::get_size(Side2D side) const {
    return axes1(side) ? m_size1 : m_size2;
}

void Block::set_size(Side2D side, int N) {
    if (axes1(side)) {
        m_size1 = N;
    } else {
        m_size2 = N;
    }

    for (Side2D l: {side, opp(side)}) {
        if (!m_adjacent_blocks[l].expired()) {
            BaseNode::Ref v1 = m_base_nodes[node(l, 0)];
            BaseNode::Ref v2 = m_base_nodes[node(l, 1)];
            if (m_adjacent_blocks[l].lock()->size(v1, v2) != N) {
                m_adjacent_blocks[l].lock()->set_size(v1, v2, N);
            }
        }
    }
}

void Block::set_size(BaseNode::Ref v1, BaseNode::Ref v2, int N) {
    set_size(get_side(v1, v2), N);
}

BsVertex::Ptr &Block::corner_vertex(int v_idx) {
    switch (v_idx) {
        case 0: return m_vertices( 0,  0);
        case 1: return m_vertices(-1,  0);
        case 2: return m_vertices(-1, -1);
        case 3: return m_vertices( 0, -1);
        default: throw std::runtime_error("Invalid vertex index");
    }
}

BsVertex::Ptr &Block::boundary_vertex(Side2D side, int idx) {
    switch (side) {
        case Side2D::B: return m_vertices(idx, 0);
        case Side2D::R: return m_vertices(-1, idx);
        case Side2D::T: return m_vertices(m_size1 - idx, -1);
        case Side2D::L: return m_vertices(0, m_size2 - idx);
        default: throw std::runtime_error("Invalid vertex index");
    }
}

BsVertex::Ptr &Block::preboundary_vertex(Side2D side, int idx) {
    switch (side) {
        case Side2D::B: return m_vertices(idx, 1);
        case Side2D::R: return m_vertices(-2, idx);
        case Side2D::T: return m_vertices(m_size1 - idx, -2);
        case Side2D::L: return m_vertices(1, m_size2 - idx);
        default: throw std::runtime_error("Invalid vertex index");
    }
}

Side2D Block::neib_face(Side2D side) const {
    int f_idx = to_face_idx(side);
    return idx_to_side((f_idx - m_rotations[side] + 4) % 4);
}

constexpr int rotation(Side2D side1, Side2D side2) {
    int f1 = to_face_idx(side1);
    int f2 = to_face_idx(side2);
    return (f1 - f2 + 4) % 4;
}

void Block::link(Block::Ref B1, Block::Ref B2) {
    if (!B1 || !B2) {
        throw std::runtime_error("Invalid block");
    }
    if (B1 == B2) {
        throw std::runtime_error("Can't link same blocks");
    }

    // Ищем просто полным перебором
    for (Side2D side1: Side2D::items()) {
        BaseNode::Ptr a1 = B1->base_node(node(side1, 0));
        BaseNode::Ptr b1 = B1->base_node(node(side1, 1));
        if (a1.get() > b1.get()) std::swap(a1, b1);

        for (Side2D side2: Side2D::items()) {
            BaseNode::Ptr a2 = B2->base_node(node(side2, 0));
            BaseNode::Ptr b2 = B2->base_node(node(side2, 1));
            if (a2.get() > b2.get()) std::swap(a2, b2);

            if (a1 == a2 && b1 == b2) {
                if (B1->m_boundaries[side1]) {
                    throw std::runtime_error("Try link through boundary #1");
                }
                B1->m_adjacent_blocks[side1] = B2;
                B1->m_rotations[side1] = rotation(side1, side2);

                if (B2->m_boundaries[side2]) {
                    throw std::runtime_error("Try link through boundary #2");
                }
                B2->m_adjacent_blocks[side2] = B1;
                B2->m_rotations[side2] = rotation(side2, side1);
            }
        }
    }
}

void Block::create_vertices() {
    using geom::SqQuad;

    m_vertices.resize(m_size1 + 1, m_size2 + 1, nullptr);

    Vector3d v0 = base_node(0)->pos();
    Vector3d v1 = base_node(1)->pos();
    Vector3d v2 = base_node(2)->pos();
    Vector3d v3 = base_node(3)->pos();
    Vector3d vL = center(Side2D::L);
    Vector3d vR = center(Side2D::R);
    Vector3d vB = center(Side2D::B);
    Vector3d vT = center(Side2D::T);

    Vector3d vC = (vL + vR + vB + vT) / 4.0;
    SqQuad quad = {
            v0, vB, v1,
            vL, vC, vR,
            v2, vT, v3
    };

    // На границах блока (низ и верх)
    for (int i = 0; i <= m_size1; ++i) {
        double x = (2.0 * i - m_size1) / m_size1;
        m_vertices(i,  0) = BsVertex::create(quad(x, -1.0));
        m_vertices(i, -1) = BsVertex::create(quad(x, +1.0));
    }

    // На границах блока (слева и справа)
    for (int j = 0; j <= m_size2; ++j) {
        double y = (2.0 * j - m_size2) / m_size2;
        m_vertices( 0, j) = BsVertex::create(quad(-1.0, y));
        m_vertices(-1, j) = BsVertex::create(quad(+1.0, y));
    }

    // Сглаживание границы
    for (int k = 0; k < 10; ++k) {
        for (Side2D side: Side2D::items()) {
            if (!boundary(side)) { continue; }

            for (int idx = 1; idx < get_size(side); ++idx) {
                Vector3d v = (boundary_vertex(side, idx - 1)->v1 +
                              boundary_vertex(side, idx + 1)->v1) / 2.0;
                v = boundary(side)->projection(v);
                boundary_vertex(side, idx)->v2 = v;
            }
        }

        for (Side2D side: Side2D::items()) {
            if (!boundary(side)) { continue; }

            for (int idx = 1; idx < get_size(side); ++idx) {
                boundary_vertex(side, idx)->v1 = boundary_vertex(side, idx)->v2;
            }
        }
    }

    /*
    // Генерация внутренних вершин, билинейное отображение
    for (int i = 1; i < m_size1; ++i) {
        for (int j = 1; j < m_size2; ++j) {
            double x = (2.0 * i - m_size1) / m_size1;
            double y = (2.0 * j - m_size2) / m_size2;
            m_vertices(i, j) = BsVertex::create(quad(x, y));
        }
    }
    */

    // Генерация внутренних вершин, Coons patch.
    for (int i = 1; i < m_size1; ++i) {
        for (int j = 1; j < m_size2; ++j) {
            double t = double(i) / m_size1;
            double s = double(j) / m_size2;

            Vector3d Lc = (1 - t) * m_vertices(0, j)->v1 + t * m_vertices(-1, j)->v1;
            Vector3d Ld = (1 - s) * m_vertices(i, 0)->v1 + s * m_vertices(i, -1)->v1;
            Vector3d B = (1 - t) * (1 - s) * m_vertices(0, 0)->v1 + (1 - t) * s * m_vertices(0, -1)->v1 +
                         t * (1 - s) * m_vertices(-1, 0)->v1 + t * s * m_vertices(-1, -1)->v1;
            Vector3d v = Lc + Ld - B;
            m_vertices(i, j) = BsVertex::create(v.x(), v.y());
        }
    }

    // Связать вершины с соседом
    for (Side2D side: Side2D::items()) {
        auto neib = m_adjacent_blocks[side];

        if (neib.expired()) { continue; }

        // Сосед есть и уже с вершинами
        if (neib.lock()->m_vertices.empty()) {
            continue;
        }

        Side2D twin = neib_face(side);

        int N = get_size(side);
        for (int k = 0; k <= N; ++k) {
            boundary_vertex(side, k) = neib.lock()->boundary_vertex(twin, N - k);
        }
    }
}

void Block::link_vertices() {
    // Базисные вершины
    for (int v_idx = 0; v_idx < 4; ++v_idx) {
        auto bv = m_base_nodes[v_idx];
        auto corner = corner_vertex(v_idx);
        corner->fix();

        /// Фиксированная точка
        if (bv->is_fixed()) { continue; }

        // Угловая вершина области, должна быть неподвижной
        if (bv->n_adjacent_blocks() < 2) { continue; }

        // Угловая вершина области по другому критерию
        auto [side1, side2] = node_adjacent_sides(v_idx);
        if (m_boundaries[side1] && m_boundaries[side2]) {
            continue;
        }

        std::vector<BsVertex::Ptr> adj_vertices;

        auto& adj_nodes = bv->adjacent_nodes();
        auto& adj_blocks = bv->adjacent_blocks();

        if (index == 2 && v_idx == 0) {

            std::cout << "Hello\n";
        }

        // Вершина на границе
        if (bv->is_boundary()) {
            if (adj_blocks.size() + 1 != adj_nodes.size()) {
                throw std::runtime_error("Block::link_vertices: n_blocks + 1 != n_nodes for boundary node");
            }
            for (int i = 0; i < adj_blocks.size(); ++i) {
                BaseNode::Ptr v1 = adj_nodes[i].lock();
                Block::Ptr block = adj_blocks[i].lock();
                if (!v1 || !block) { throw std::runtime_error("Block::link_vertices: invalid block"); }
                auto side = block->get_side(bv, v1);
                adj_vertices.emplace_back(block->boundary_vertex(side, 1));

                if (block->boundary(side)) {
                    corner->add_boundary(block->boundary(side).get());
                }
            }
            BaseNode::Ptr v1 = adj_nodes.back().lock();
            Block::Ptr block = adj_blocks.back().lock();
            if (!v1 || !block) { throw std::runtime_error("Block::link_vertices: invalid block"); }
            auto side = block->get_side(bv, v1);
            adj_vertices.emplace_back(block->boundary_vertex(side, block->get_size(side) - 1));
            if (block->boundary(side)) {
                corner->add_boundary(block->boundary(side).get());
            }
        } else {
            if (adj_nodes.size() != adj_blocks.size()) {
                throw std::runtime_error("Block::link_vertices: n_blocks != n_nodes for inner node");
            }
            for (int i = 0; i < adj_blocks.size(); ++i) {
                BaseNode::Ptr v1 = adj_nodes[i].lock();
                Block::Ptr block = adj_blocks[i].lock();
                if (!v1 || !block) { throw std::runtime_error("Block::link_vertices: invalid block"); }
                auto side = block->get_side(bv, v1);
                adj_vertices.emplace_back(block->boundary_vertex(side, 1));
            }
        }
        corner->set_adjacent_vertices(adj_vertices);
    }

    // Вершины на границах блока
    for (Side2D side: Side2D::items()) {
        if (!boundary(side) && !adjacent_block(side)) {
            throw std::runtime_error("BFace is not boundary and has no neighbor");
        }

        // Число ячеек по грани
        if (boundary(side)) {
            for (int idx = 1; idx < get_size(side); ++idx) {
                boundary_vertex(side, idx)->set_adjacent_vertices(
                    {
                        boundary_vertex(side, idx + 1),
                        preboundary_vertex(side, idx),
                        boundary_vertex(side, idx - 1)
                    }
                );
                boundary_vertex(side, idx)->add_boundary(
                    m_boundaries[side].get()
                );
            }
        } else {
            int N = get_size(side);
            auto neib = m_adjacent_blocks[side].lock();
            Side2D twin = neib_face(side);

            for (int idx = 1; idx < N; ++idx) {
                auto K0 = boundary_vertex(side, idx);
                auto K1 = boundary_vertex(side, idx - 1);
                auto K2 = boundary_vertex(side, idx + 1);
                auto K3 = preboundary_vertex(side, idx);
                auto K4 = neib->preboundary_vertex(twin, N - idx);

                boundary_vertex(side, idx)->set_adjacent_vertices(
                    {
                        boundary_vertex(side, idx - 1),
                        boundary_vertex(side, idx + 1),
                        preboundary_vertex(side, idx),
                        neib->preboundary_vertex(twin, N - idx)
                    });
            }
        }
    }

    // Внутренние вершины
    for (int i = 1; i < m_size1; ++i) {
        for (int j = 1; j < m_size2; ++j) {
            m_vertices(i, j)->set_adjacent_vertices(
                {
                    m_vertices(i + 1, j),
                    m_vertices(i, j + 1),
                    m_vertices(i - 1, j),
                    m_vertices(i, j - 1)
                }
            );
        }
    }
}

double Block::calc_modulus() const {
    int Nx = m_size1;
    int Ny = m_size2;

    double sum_x = 0.0;
    double sum_y = 0.0;
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            Vector3d v1 = m_vertices(i, j)->v1;
            Vector3d v2 = m_vertices(i+1, j)->v1;
            Vector3d v3 = m_vertices(i+1, j+1)->v1;
            Vector3d v4 = m_vertices(i, j+1)->v1;

            // Triangle Method
            Triangle T1 = {v4, v1, v2};
            Triangle T2 = {v1, v2, v3};
            Triangle T3 = {v2, v3, v4};
            Triangle T4 = {v3, v4, v1};

            double A1 = T1.area();
            double A2 = T2.area();
            double A3 = T3.area();
            double A4 = T4.area();

            double int_x = (v2 - v3).squaredNorm() * (1/A2 + 1/A3) + (v1 - v4).squaredNorm() * (1/A1 + 1/A4);
            double int_y = (v1 - v2).squaredNorm() * (1/A1 + 1/A2) + (v3 - v4).squaredNorm() * (1/A3 + 1/A4);
            int_x /= (8 * Nx * Nx);
            int_y /= (8 * Ny * Ny);

            sum_x += int_x;
            sum_y += int_y;
        }
    }

    std::cout << "M: " << sum_y << "; 1/M: " << sum_x << "; 1: " << sum_x * sum_y << std::endl;
    std::cout << "M avg: " << std::sqrt(sum_y / sum_x) << std::endl;
    std::cout << "M avg: " << 0.5 * (sum_y + 1.0 / sum_x) << std::endl;
    return std::sqrt(sum_y / sum_x);
}

double Block::smooth() {
    double err = 0.0;

    for (int i = 0; i <= size1(); ++i) {
        for (int j = 0; j <= size2(); ++j) {
            auto vertex = m_vertices(i, j);
            auto &adjacent = vertex->adjacent_vertices();

            // Фиксированная вершина
            if (adjacent.empty()) {
                vertex->v2 =vertex->v1;
                continue;
            }

            // Угловая вершина
            if (!vertex->inner() && vertex->corner()) {
                vertex->v2 = vertex->v1;
                continue;
            }

            // Вершина на границе
            if (!vertex->inner()) {
                Curve *boundary = vertex->boundary();

                double K =  (double(size2()) / double(size1()));

                // Какому направлению блока принадлежит
                int dir = (m_boundaries[Side2D::B].get() == boundary ||
                           m_boundaries[Side2D::T].get() == boundary) ? 0 : 1;

                double sum = 0.0;
                double w1 = 0.5;
                double w2 = 1.0;

                if (dir == 0) {
                    w1 /= (m_modulus * K);
                    w2 *= (m_modulus * K);
                }
                else {
                    w1 *= (m_modulus * K);
                    w2 /= (m_modulus * K);
                }

                Vector3d next = Vector3d::Zero();
                next += w1 * adjacent[0]->v1;
                sum += w1;
                for (int n = 1; n < int(adjacent.size()) - 1; ++n) {
                    next += w2 * boundary->projection(adjacent[n]->v1);
                    sum += w2;
                }
                sum += w1;
                next += w1 * adjacent.back()->v1;
                next /= sum;

                vertex->v2 = vertex->boundary()->projection(next);
                continue;
            }

            // Внутренняя вершина
            // А может по отдельности все зафигачить???
            // Есть 4 случая:
            //   Регулярная внутренняя
            //   Регулярная на границе двух блоков
            //   Регулярная на границе 4ех блоков
            //   Сингулярная на границе нескольких блоков

            // Проигнорируем сингулярные (пока)
            if (adjacent.size() != 4) {
                vertex->v2 = vertex->v1;
                throw std::runtime_error("NO WAY #0");
                continue;
            }

            double ML, MR, MB, MT;
            double NL, NR, NB, NT;

            // Классический вариант внутри блоков
            MB = MT = modulus();
            ML = MR = 1.0 / modulus();
            NB = NT = ratio();
            NL = NR = 1.0 / ratio();

            BsVertex* vR = adjacent[0];
            BsVertex* vT = adjacent[1];
            BsVertex* vL = adjacent[2];
            BsVertex* vB = adjacent[3];

            // Регулярная внутри блока
            if (0 < i && i < size1() && 0 < j && j < size2()) {
                vL = m_vertices(i - 1, j).get();
                vR = m_vertices(i + 1, j).get();
                vB = m_vertices(i, j - 1).get();
                vT = m_vertices(i, j + 1).get();
            }
            // На нижней границе блока
            else if (0 < i && i < size1() && j == 0) {
                vL = m_vertices(i - 1, j).get();
                vR = m_vertices(i + 1, j).get();
                vT = m_vertices(i, j + 1).get();

                auto neib = m_adjacent_blocks[Side2D::B].lock();
                vB = neib->vertex(i, -2).get();
                MB = neib->modulus();
                NB = neib->ratio();

                ML = MR = 1.0 / std::sqrt(modulus() * neib->modulus());
                NL = NR = 1.0 / std::sqrt(ratio() * neib->ratio());
            }
            // На верхней границе блока
            else if (0 < i && i < size1() && j == size2()) {
                vL = m_vertices(i - 1, j).get();
                vR = m_vertices(i + 1, j).get();
                vB = m_vertices(i, j - 1).get();

                auto neib = m_adjacent_blocks[Side2D::T].lock();
                vT = neib->vertex(i, 1).get();
                MT = neib->modulus();
                NT = neib->ratio();

                ML = MR = 1.0 / std::sqrt(modulus() * neib->modulus());
                NL = NR = 1.0 / std::sqrt(ratio() * neib->ratio());
            }
            // На левой границе блока
            else if (0 < j && j < size2() && i == 0) {
                vR = m_vertices(i + 1, j).get();
                vB = m_vertices(i, j - 1).get();
                vT = m_vertices(i, j + 1).get();

                auto neib = m_adjacent_blocks[Side2D::L].lock();
                vL = neib->vertex(-2, j).get();
                ML = 1.0 / neib->modulus();
                NL = 1.0 / neib->ratio();

                MB = MT = std::sqrt(modulus() * neib->modulus());
                NB = NT = std::sqrt(ratio() * neib->ratio());
            }
            // На правой границе блока
            else if (0 < j && j < size2() && i == size1()) {
                vL = m_vertices(i - 1, j).get();
                vB = m_vertices(i, j - 1).get();
                vT = m_vertices(i, j + 1).get();

                auto neib = m_adjacent_blocks[Side2D::R].lock();
                vR = neib->vertex(1, j).get();
                MR = 1.0 / neib->modulus();
                NR = 1.0 / neib->ratio();

                MB = MT = std::sqrt(modulus() * neib->modulus());
                NB = NT = std::sqrt(ratio() * neib->ratio());
            }
            else {
                // Регулярная угловая точка
                // У меня сейчас одна точка, рассмотрим её в лоб, для конкретного блока
                if (i != 0 || j != 0) { continue; }

                auto rt = shared_from_this();
                auto lt = m_adjacent_blocks[Side2D::L].lock();
                auto rb = m_adjacent_blocks[Side2D::B].lock();
                auto lb = lt->m_adjacent_blocks[Side2D::B].lock();

                vR = m_vertices(1, 0).get();
                vT = m_vertices(0, 1).get();
                vL = lt->vertex(-2, 0).get();
                vB = rb->vertex(0, -2).get();

                MB = std::sqrt(lb->modulus() * rb->modulus());
                MT = std::sqrt(lt->modulus() * rt->modulus());
                ML = 1.0 / std::sqrt(lb->modulus() * lt->modulus());
                MR = 1.0 / std::sqrt(rb->modulus() * rt->modulus());
                
                NB = std::sqrt(lb->ratio() * rb->ratio());
                NT = std::sqrt(lt->ratio() * rt->ratio());
                NL = 1.0 / std::sqrt(lb->ratio() * lt->ratio());
                NR = 1.0 / std::sqrt(rb->ratio() * rt->ratio());
            }

            double KL = ML * NL;
            double KR = MR * NR;
            double KB = MB * NB;
            double KT = MT * NT;

            Vector3d next = Vector3d::Zero();
            next = KL * vL->v1 + KR * vR->v1 + KB * vB->v1 + KT * vT->v1;
            next /= (KL + KR + KB + KT);

            // for (auto neib: adjacent) {
            //     next += neib->v1;
            // }
            // next /= vertex->n_adjacent();

            vertex->v2 = next;

            double L = 0.0;
            for (const auto adj: adjacent) {
                L = std::max(L, (vertex->v1 - adj->v1).norm());
            }

            err = std::max(err, (vertex->v1 - vertex->v2).norm() / L);
        }
    }

    return err;
}

void Block::update() {
    for (int i = 0; i < m_vertices.size1(); ++i) {
        for (int j = 0; j < m_vertices.size2(); ++j) {
            m_vertices(i, j)->v1 = m_vertices(i, j)->v2;
        }
    }
}

} // namespace zephyr::geom::generator
