#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <format>
#include <bits/fs_fwd.h>

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

constexpr int node(Side2D side, int idx) {
    switch (side) {
        case Side2D::L: return std::array{0, 2}[idx];
        case Side2D::R: return std::array{1, 3}[idx];
        case Side2D::B: return std::array{0, 1}[idx];
        case Side2D::T: return std::array{2, 3}[idx];
        default: throw std::runtime_error("Invalid side");
    }
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
    switch (base_node_index(v)) {
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
    return node_adjacent_sides(base_node_index(v));
}

std::tuple<Side2D, Side2D> Block::adjacent_sides(BaseNode::Ref v) const {
    return adjacent_sides(v.get());
}

Block::Ptr Block::adjacent_block(Side2D side) const {
    return m_adjacent_blocks[side].lock();
}

int Block::base_node_index(const BaseNode* v) const {
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
    int idx_v1 = base_node_index(v1);
    int idx_v2 = base_node_index(v2);
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

int Block::get_axis(BaseNode::Ref v1, BaseNode::Ref v2) const {
    auto side = get_side(v1, v2);
    return to_axis(side);
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

void Block::set_boundary(Side2D side, Curve::Ref curve) {
    m_boundaries[side] = curve;
    m_adjacent_blocks[side].reset();
}

void Block::set_boundary(BaseNode::Ref v1, BaseNode::Ref v2, Curve::Ref curve) {
    set_boundary(get_side(v1, v2), curve);
}

void Block::set_size(int axis, int N) {
    m_sizes[axis] = N;
    for (Side2D l: sides(axis)) {
        if (!m_adjacent_blocks[l].expired()) {
            BaseNode::Ref v1 = m_base_nodes[node(l, 0)];
            BaseNode::Ref v2 = m_base_nodes[node(l, 1)];
            if (m_adjacent_blocks[l].lock()->size(v1, v2) != N) {
                m_adjacent_blocks[l].lock()->set_size(v1, v2, N);
            }
        }
    }
}

void Block::set_size(Side2D side, int N) {
    set_size(to_axis(side), N);
}

void Block::set_size(BaseNode::Ref v1, BaseNode::Ref v2, int N) {
    set_size(get_axis(v1, v2), N);
}

BsVertex::Ptr& Block::corner_vertex(int v_idx) {
    switch (v_idx) {
        case 0: return vertex( 0,  0);
        case 1: return vertex(-1,  0);
        case 2: return vertex( 0, -1);
        case 3: return vertex(-1, -1);
        default: throw std::runtime_error("Invalid vertex index");
    }
}

BsVertex::Ptr& Block::boundary_vertex(Side2D side, int idx) {
    switch (side) {
        case Side2D::B: return vertex(idx, 0);
        case Side2D::R: return vertex(-1, idx);
        case Side2D::T: return vertex(size1() - idx, -1);
        case Side2D::L: return vertex(0, size2() - idx);
        default: throw std::runtime_error("Invalid vertex index");
    }
}

BsVertex::Ptr& Block::preboundary_vertex(Side2D side, int idx) {
    switch (side) {
        case Side2D::B: return vertex(idx, 1);
        case Side2D::R: return vertex(-2, idx);
        case Side2D::T: return vertex(size1() - idx, -2);
        case Side2D::L: return vertex(1, size2() - idx);
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
                if (B1->boundary(side1)) {
                    throw std::runtime_error("Try link through boundary #1");
                }
                B1->m_adjacent_blocks[side1] = B2;
                B1->m_rotations[side1] = rotation(side1, side2);

                if (B2->boundary(side2)) {
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

    m_vertices.resize({size1() + 1, size2() + 1}, nullptr);

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
    for (int i = 0; i <= size1(); ++i) {
        double x = (2.0 * i - size1()) / size1();
        vertex(i,  0) = BsVertex::create(quad(x, -1.0));
        vertex(i, -1) = BsVertex::create(quad(x, +1.0));
    }

    // На границах блока (слева и справа)
    for (int j = 0; j <= size2(); ++j) {
        double y = (2.0 * j - size2()) / size2();
        vertex( 0, j) = BsVertex::create(quad(-1.0, y));
        vertex(-1, j) = BsVertex::create(quad(+1.0, y));
    }

    // Сглаживание границы
    for (int k = 0; k < 10; ++k) {
        for (Side2D side: Side2D::items()) {
            if (!boundary(side)) { continue; }

            for (int idx = 1; idx < size(side); ++idx) {
                Vector3d v = (boundary_vertex(side, idx - 1)->v1 +
                              boundary_vertex(side, idx + 1)->v1) / 2.0;
                v = boundary(side)->projection(v);
                boundary_vertex(side, idx)->v2 = v;
            }
        }

        for (Side2D side: Side2D::items()) {
            if (!boundary(side)) { continue; }

            for (int idx = 1; idx < size(side); ++idx) {
                boundary_vertex(side, idx)->v1 = boundary_vertex(side, idx)->v2;
            }
        }
    }

    /*
    // Генерация внутренних вершин, биквадратичное отображение
    for (int i = 1; i < size1(); ++i) {
        for (int j = 1; j < size2(); ++j) {
            double x = (2.0 * i - size1()) / size1();
            double y = (2.0 * j - size2()) / size2();
            vertex(i, j) = BsVertex::create(quad(x, y));
        }
    }
    */

    // Генерация внутренних вершин, Coons patch.
    for (int i = 1; i < size1(); ++i) {
        for (int j = 1; j < size2(); ++j) {
            double t = double(i) / size1();
            double s = double(j) / size2();

            Vector3d Lc = (1 - t) * vertex(0, j)->v1 + t * vertex(-1, j)->v1;
            Vector3d Ld = (1 - s) * vertex(i, 0)->v1 + s * vertex(i, -1)->v1;
            Vector3d B = (1 - t) * (1 - s) * vertex(0, 0)->v1 + (1 - t) * s * vertex(0, -1)->v1 +
                         t * (1 - s) * vertex(-1, 0)->v1 + t * s * vertex(-1, -1)->v1;
            Vector3d v = Lc + Ld - B;
            vertex(i, j) = BsVertex::create(v.x(), v.y());
        }
    }

    // Связать вершины с соседом
    for (Side2D side: Side2D::items()) {
        auto neib = adjacent_block(side);

        if (!neib) { continue; }

        // Сосед есть, но без вершин
        if (neib->vertices().empty()) {
            continue;
        }

        Side2D twin = neib_face(side);

        int N = size(side);
        for (int k = 0; k <= N; ++k) {
            boundary_vertex(side, k) = neib->boundary_vertex(twin, N - k);
        }
    }
}

void Block::link_vertices() {
    // Очистить все списки
    for (int i = 0; i <= size1(); ++i) {
        for (int j = 0; j <= size2(); ++j) {
            vertex(i, j)->clear();
        }
    }

    // Внутренние вершины блоков
    for (int i = 1; i < size1(); ++i) {
        for (int j = 1; j < size2(); ++j) {
            // Справа, сверху, слева, снизу
            std::vector<BsVertex::Edge> edges {
                {.neib = vertex(i + 1, j).get(), .ratio1 = ratio_ptr(0), .ratio2 = ratio_ptr(0)},
                {.neib = vertex(i, j + 1).get(), .ratio1 = ratio_ptr(1), .ratio2 = ratio_ptr(1)},
                {.neib = vertex(i - 1, j).get(), .ratio1 = ratio_ptr(0), .ratio2 = ratio_ptr(0)},
                {.neib = vertex(i, j - 1).get(), .ratio1 = ratio_ptr(1), .ratio2 = ratio_ptr(1)}
            };
            vertex(i, j)->set_edges(edges);
        }
    }

    // Вершины на границах блока (без угловых)
    for (Side2D side: Side2D::items()) {
        if (!boundary(side) && !adjacent_block(side)) {
            throw std::runtime_error("BFace is not boundary and has no neighbor");
        }

        int axis = to_axis(side);
        int perp_axis = (axis + 1) % 2;

        if (boundary(side)) {
            for (int idx = 1; idx < size(side); ++idx) {
                // Добавить границу
                boundary_vertex(side, idx)->add_boundary(boundary(side).get());

                // Добавить связи (три штуки)
                std::vector<BsVertex::Edge> edges;
                edges.push_back({
                    .neib   = boundary_vertex(side, idx + 1).get(),
                    .ratio1 = ratio_ptr(axis),
                    .ratio2 = nullptr
                });
                edges.push_back({
                    .neib   = preboundary_vertex(side, idx).get(),
                    .ratio1 = ratio_ptr(perp_axis),
                    .ratio2 = ratio_ptr(perp_axis)
                });
                edges.push_back({
                    .neib   = boundary_vertex(side, idx - 1).get(),
                    .ratio1 = ratio_ptr(axis),
                    .ratio2 = nullptr
                });
                boundary_vertex(side, idx)->set_edges(edges);
            }
        } else {
            auto neib = adjacent_block(side);
            if (!neib) {
                throw std::runtime_error("Face is not boundary and has no neighbor");
            }

            int N = size(side);
            Side2D twin = neib_face(side);

            int neib_axis = to_axis(twin);
            int prep_neib_axis = (neib_axis + 1) % 2;

            for (int idx = 1; idx < N; ++idx) {
                std::vector<BsVertex::Edge> edges;
                edges.push_back({
                    .neib   = boundary_vertex(side, idx + 1).get(),
                    .ratio1 = ratio_ptr(axis),
                    .ratio2 = neib->ratio_ptr(neib_axis)
                });
                edges.push_back({
                    .neib   = preboundary_vertex(side, idx).get(),
                    .ratio1 = ratio_ptr(perp_axis),
                    .ratio2 = ratio_ptr(perp_axis)
                });
                edges.push_back({
                    .neib   = boundary_vertex(side, idx - 1).get(),
                    .ratio1 = ratio_ptr(axis),
                    .ratio2 = neib->ratio_ptr(neib_axis)
                });
                edges.push_back({
                    .neib   = neib->preboundary_vertex(twin, N - idx).get(),
                    .ratio1 = neib->ratio_ptr(prep_neib_axis),
                    .ratio2 = neib->ratio_ptr(prep_neib_axis)
                });
                boundary_vertex(side, idx)->set_edges(edges);
            }
        }
    }

    // Угловые вершины блока (здесь могут быть сингулярности)
    for (int v_idx = 0; v_idx < 4; ++v_idx) {
        auto node = base_node(v_idx);
        corner_vertex(v_idx)->clear();

        /// Фиксированная точка
        if (node->is_fixed()) { continue; }

        // Угловая вершина области, должна быть неподвижной
        if (node->n_adjacent_blocks() < 2) { continue; }

        // Угловая вершина области по другому критерию
        if (auto [side1, side2] = node_adjacent_sides(v_idx);
            m_boundaries[side1] && m_boundaries[side2]) {
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

            auto corner = corner_vertex(v_idx);
            std::vector<BsVertex::Edge> edges;

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
                edges.push_back({
                    .neib = block->boundary_vertex(side, 1).get(),
                    .ratio1 = block->ratio_ptr(to_axis(side)),
                    .ratio2 = nullptr
                });
            }
            // Внутренние вершины
            for (int i = 1; i < adjacent_blocks.size(); ++i) {
                BaseNode::Ptr node2 = adjacent_nodes[i].lock();
                if (!node2) { throw std::runtime_error("Block::link_vertices: invalid node"); }
                Block::Ptr block = adjacent_blocks[i].lock();
                auto prev_block = adjacent_blocks[i - 1].lock();
                if (!block || !prev_block) { throw std::runtime_error("Block::link_vertices: invalid block"); }
                Side2D side = block->get_side(node, node2);
                Side2D prev_side = prev_block->get_side(node, node2);
                edges.push_back({
                    .neib   = block->boundary_vertex(side, 1).get(),
                    .ratio1 = block->ratio_ptr(side),
                    .ratio2 = prev_block->ratio_ptr(prev_side),
                });
                if (block->boundary(side)) {
                    throw std::runtime_error("Block::link_vertices: invalid block (why boundary)");
                }
            }
            // Последняя граничная вершина
            {
                BaseNode::Ptr node2 = adjacent_nodes.back().lock();
                Block::Ptr block = adjacent_blocks.back().lock();
                if (!node2 || !block) { throw std::runtime_error("Block::link_vertices: invalid block"); }
                Side2D side = block->get_side(node, node2);
                edges.push_back({
                    .neib = block->boundary_vertex(side, block->size(side) - 1).get(),
                    .ratio1 = block->ratio_ptr(side),
                    .ratio2 = nullptr
                });
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

            auto corner = corner_vertex(v_idx);
            std::vector<BsVertex::Edge> edges;

            for (int i = 0; i < adjacent_blocks.size(); ++i) {
                BaseNode::Ptr node2 = adjacent_nodes[i].lock();
                if (!node2) { throw std::runtime_error("Block::link_vertices: invalid node"); }
                Block::Ptr block = adjacent_blocks[i].lock();
                int j = (i + adjacent_blocks.size() - 1) % adjacent_blocks.size();
                auto prev_block = adjacent_blocks[j].lock();
                if (!block || !prev_block) { throw std::runtime_error("Block::link_vertices: invalid block"); }
                Side2D side = block->get_side(node, node2);
                Side2D prev_side = prev_block->get_side(node, node2);
                edges.push_back({
                    .neib   = block->boundary_vertex(side, 1).get(),
                    .ratio1 = block->ratio_ptr(side),
                    .ratio2 = prev_block->ratio_ptr(prev_side),
                });
                if (block->boundary(side)) {
                    throw std::runtime_error("Block::link_vertices: invalid block (why boundary)");
                }
            }
            corner->set_edges(edges);
        }
    }
}

double Block::calc_modulus() const {
    int Nx = size1();
    int Ny = size2();

    double sum_x = 0.0;
    double sum_y = 0.0;
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            Vector3d v1 = vertex(i, j)->v1;
            Vector3d v2 = vertex(i+1, j)->v1;
            Vector3d v3 = vertex(i+1, j+1)->v1;
            Vector3d v4 = vertex(i, j+1)->v1;

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

    //std::cout << "M: " << sum_y << "; 1/M: " << sum_x << "; 1: " << sum_x * sum_y << std::endl;
    //std::cout << "M avg: " << std::sqrt(sum_y / sum_x) << std::endl;
    //std::cout << "M avg: " << 0.5 * (sum_y + 1.0 / sum_x) << std::endl;
    return std::sqrt(sum_y / sum_x);
}

void Block::update_modulus() {
    double K = calc_modulus();
    m_modulus = K;
    m_ratio[0] = m_modulus * size2() / size1();
    m_ratio[1] = 1.0 / m_ratio[0];

    std::cout << std::format("Block {}.\tK: {:.2f}; ratio1: {:.2f}; ratio2: {:.2f}\n", index, m_modulus, m_ratio[0], m_ratio[1]);
}

double Block::smooth() {
    double err = 0.0;

    for (int i = 0; i <= size1(); ++i) {
        for (int j = 0; j <= size2(); ++j) {
            auto vc = vertex(i, j);
            auto &adjacent = vc->adjacent();

            // Фиксированная вершина
            if (adjacent.empty()) {
                vc->v2 =vc->v1;
                continue;
            }

            // Угловая вершина
            if (!vc->inner() && vc->corner()) {
                vc->v2 = vc->v1;
                continue;
            }

            // Вершина на границе
            if (!vc->inner()) {
                Curve *boundary = vc->boundary();

                Vector3d next = Vector3d::Zero();
                double sum = 0.0;
                for (auto edge: vertex(i, j)->adjacent()) {
                    next += 1/edge.ratio() * edge.neib->v1;
                    sum += 1/edge.ratio();
                }
                next /= sum;
                vc->v2 = vc->boundary()->projection(next);
                continue;
            }

            // Внутренняя вершина
            // А может по отдельности все зафигачить???
            // Есть 4 случая:
            //   Регулярная внутренняя
            //   Регулярная на границе двух блоков
            //   Регулярная на границе 4ех блоков
            //   Сингулярная на границе нескольких блоков


            Vector3d next = Vector3d::Zero();
            double sum = 0.0;
            for (auto edge: vertex(i, j)->adjacent()) {
                next += 1/edge.ratio() * edge.neib->v1;
                sum += 1/edge.ratio();
            }
            next /= sum;
            vc->v2 = next;
            double L = 0.0;
            for (const auto adj: adjacent) {
                L = std::max(L, (vc->v1 - adj.neib->v1).norm());
            }

            err = std::max(err, (vc->v1 - vc->v2).norm() / L);

#if 0
            // Проигнорируем сингулярные (пока)
            if (adjacent.size() != 4) {
                vc->v2 = vc->v1;
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

            BsVertex* vR = adjacent[0].neib;
            BsVertex* vT = adjacent[1].neib;
            BsVertex* vL = adjacent[2].neib;
            BsVertex* vB = adjacent[3].neib;

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
            double sum = 0.0;
            for (auto edge: vertex(i, j)->adjacent()) {
                next += 1/edge.ratio() * edge.neib->v1;
                sum += 1/edge.ratio();
            }
            next /= sum;
            //next = KL * vL->v1 + KR * vR->v1 + KB * vB->v1 + KT * vT->v1;
           // next /= (KL + KR + KB + KT);

            // for (auto neib: adjacent) {
            //     next += neib->v1;
            // }
            // next /= vertex->n_adjacent();

            vc->v2 = next;

            double L = 0.0;
            for (const auto adj: adjacent) {
                L = std::max(L, (vc->v1 - adj.neib->v1).norm());
            }

            err = std::max(err, (vc->v1 - vc->v2).norm() / L);
#endif
        }
    }

    return err;
}

void Block::update() {
    z_assert(m_vertices.size1() == size1() + 1, "Bad size1");
    z_assert(m_vertices.size2() == size2() + 1, "Bad size2");
    for (int i = 0; i <= size1(); ++i) {
        for (int j = 0; j <= size2(); ++j) {
            vertex(i, j)->v1 = vertex(i, j)->v2;
        }
    }
}

} // namespace zephyr::geom::generator
