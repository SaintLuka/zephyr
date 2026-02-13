#include <algorithm>
#include <array>
#include <map>
#include <memory>
#include <set>
#include <format>

#include <zephyr/geom/primitives/cube.h>

#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>

#include <zephyr/geom/intersection.h>
#include <zephyr/geom/primitives/triangle.h>
#include <zephyr/math/funcs.h>

#include "zephyr/geom/curves/cubic_spline.h"

namespace zephyr::geom::generator {

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

// Выбрать ось блока, вдоль которой лежит сторона
constexpr Axis to_axis(Side side) {
    return (side == Side::B || side == Side::T) ? Axis::X : Axis::Y;
}

// Стороны, которые лежат вдоль выбранной оси
constexpr std::array<Side, 2> sides_by_axis(Axis axis) {
    return axis == Axis::X? std::array{Side::B, Side::T} : std::array{Side::L, Side::R};
}

// Переход к индексации сторон против часовой стрелки
constexpr int side_to_idx(Side side) {
    switch (side) {
        case Side::BOTTOM: return 0;
        case Side::RIGHT:  return 1;
        case Side::TOP:    return 2;
        case Side::LEFT:   return 3;
        default: throw std::invalid_argument("Invalid side");
    }
}

// Переход от индексации сторон против часовой стрелки
constexpr Side idx_to_side(int idx) {
    switch (idx) {
        case 0: return Side::BOTTOM;
        case 1: return Side::RIGHT;
        case 2: return Side::TOP;
        case 3: return Side::LEFT;
        default: throw std::invalid_argument("Invalid side");
    }
}

// Индекс базисной вершины на стороне
constexpr int node_idx(Side side, int idx) {
    switch (side) {
        case Side::LEFT:   return std::array{0, 2}[idx];
        case Side::RIGHT:  return std::array{1, 3}[idx];
        case Side::BOTTOM: return std::array{0, 1}[idx];
        case Side::TOP:    return std::array{2, 3}[idx];
        default: throw std::runtime_error("Invalid side");
    }
}

// Поворот одного блока относительно другого при смежных сторонах (side1, side2)
constexpr int rotation(Side side1, Side side2) {
    int f1 = side_to_idx(side1);
    int f2 = side_to_idx(side2);
    return (f1 - f2 + 4) % 4;
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

Block::Block(const std::array<BaseNode::Ptr, 4>& vertices) {
    // Центр блока
    Vector3d vc = Vector3d::Zero();
    for (const auto& v: vertices) {
        vc += v->pos();
    }
    vc *= 0.25;

    // Пары (угол, вершина)
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
}

Block::Ptr Block::create(const std::array<BaseNode::Ptr, 4>& vertices) {
    return std::make_shared<Block>(vertices);
}

void Block::set_index(int i) {
    if (i < 0) {
        throw std::invalid_argument("Block::set_index error: Invalid index");
    }
    m_index = i;
}

// ---------- Геометрия ------------------------------------------------------------------------------------------------

Vector3d Block::center() const {
    return 0.25 * (m_base_nodes[0]->pos() + m_base_nodes[1]->pos() + m_base_nodes[2]->pos() + m_base_nodes[3]->pos());
}

Vector3d Block::center(Side side) const {
    Vector3d v1 = base_node(node_idx(side, 0))->pos();
    Vector3d v2 = base_node(node_idx(side, 1))->pos();
    Vector3d vc = 0.5 * (v1 + v2);
    if (boundary(side)) {
        vc = boundary(side)->projection(vc);
    }
    return vc;
}

double Block::length(Side side) const {
    Vector3d v1 = base_node(node_idx(side, 0))->pos();
    Vector3d v2 = base_node(node_idx(side, 1))->pos();
    return (v2 - v1).norm();
}

// ---------- Базисные вершины -----------------------------------------------------------------------------------------

int Block::base_node_index(const BaseNode* v) const {
    for (int idx = 0; idx < 4; ++idx) {
        if (m_base_nodes[idx].get() == v) {
            return idx;
        }
    }
    throw std::runtime_error("Can't find base node");
}

std::tuple<Side, Side> Block::adjacent_sides(const BaseNode* v) const {
    if (!v) {
        throw std::runtime_error("Block::adjacent_sides() error: Invalid BaseNode");
    }
    return node_adjacent_sides(base_node_index(v));
}

std::tuple<Side, Side> Block::adjacent_sides(BaseNode::Ref v) const {
    return adjacent_sides(v.get());
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
        default: throw std::runtime_error("Invalid base node index");
    }
}

std::tuple<BaseNode::Ptr, BaseNode::Ptr> Block::base_nodes(Side side) const {
    return {base_node(node_idx(side, 0)), base_node(node_idx(side, 1))};
}

// ---------- Размеры блока ----------------------------------------- --------------------------------------------------

int Block::size(Axis axis) const {
    return m_sizes[static_cast<int>(axis)];
}

int Block::size(Side side) const {
    return m_sizes[static_cast<int>(to_axis(side))];
}

int Block::size(BaseNode::Ref v1, BaseNode::Ref v2) const {
    return m_sizes[static_cast<int>(get_axis(v1, v2))];
}

void Block::set_size(Axis axis, int N) {
    if (m_sizes[static_cast<int>(axis)] == N) return;
    m_sizes[static_cast<int>(axis)] = N;
    update_lambda();
    for (Side side: sides_by_axis(axis)) {
        auto block = adjacent_block(side);
        if (block) {
            BaseNode::Ref v1 = m_base_nodes[node_idx(side, 0)];
            BaseNode::Ref v2 = m_base_nodes[node_idx(side, 1)];
            if (block->size(v1, v2) != N) {
                block->set_size(v1, v2, N);
            }
        }
    }
}

void Block::set_size(Side side, int N) {
    set_size(to_axis(side), N);
}

void Block::set_size(BaseNode::Ref v1, BaseNode::Ref v2, int N) {
    set_size(get_axis(v1, v2), N);
}

std::string Block::sizes_info() const {
    return std::format("Block {:3}.  K: {:.2f};  sizes: ({:3}, {:3})", m_index, m_modulus, m_sizes[0], m_sizes[1]);
}

double Block::rel_size(Axis axis) const {
    return m_rel_sizes[static_cast<int>(axis)];
}

double Block::rel_size(BaseNode::Ref v1, BaseNode::Ref v2) const {
    return m_rel_sizes[static_cast<int>(get_axis(v1, v2))];
}

void Block::reset_rel_sizes() {
    m_rel_sizes[0] = NAN;
    m_rel_sizes[1] = NAN;
}

void Block::init_rel_sizes() {
    if (std::isnan(m_modulus)) {
        throw std::runtime_error("Invalid Modulus");
    }
    set_rel_size(Axis::X, m_modulus);
    set_rel_size(Axis::Y, 1.0);
}

bool Block::update_rel_sizes() {
    if (std::isnan(m_modulus)) {
        throw std::runtime_error("update_rel_sizes error: Invalid Modulus");
    }
    if (bad_number(rel_size(Axis::X)) && good_number(rel_size(Axis::Y))) {
        set_rel_size(Axis::X, rel_size(Axis::Y) * modulus());
    }
    if (bad_number(rel_size(Axis::Y)) && good_number(rel_size(Axis::X))) {
        set_rel_size(Axis::Y, rel_size(Axis::X) / modulus());
    }
    return good_number(rel_size(Axis::X)) && good_number(rel_size(Axis::Y));
}

void Block::set_rel_size(Axis axis, double N) {
    if (m_rel_sizes[static_cast<int>(axis)] == N) return;

    m_rel_sizes[static_cast<int>(axis)] = N;
    //std::cout << "  Block " << index << ". set size axis " << static_cast<int>(axis) << ": " << N << "; ";
    //std::cout << "new sizes: (" << m_rel_sizes[0] << ", " << m_rel_sizes[1] << ")\n";
    for (Side side: sides_by_axis(axis)) {
        auto block = adjacent_block(side);
        if (block) {
            BaseNode::Ref v1 = m_base_nodes[node_idx(side, 0)];
            BaseNode::Ref v2 = m_base_nodes[node_idx(side, 1)];
            if (block->rel_size(v1, v2) != N) {
                block->set_rel_size(v1, v2, N);
            }
        }
    }
}

void Block::set_rel_size(Side side, double N) {
    set_rel_size(to_axis(side), N);
}

void Block::set_rel_size(BaseNode::Ref v1, BaseNode::Ref v2, double N) {
    set_rel_size(get_axis(v1, v2), N);
}

void Block::set_rel_sizes(double Nx, double Ny) {
    set_rel_size(Axis::X, Nx);
    set_rel_size(Axis::Y, Ny);
}

// ---------- Границы блока ------------------ -------------------------------------------------------------------------

Curve::Ref Block::boundary(Side side) const {
    return m_boundaries[static_cast<int>(side)];
}

bool Block::is_boundary(const BaseNode* v) const {
    auto [side1, side2] = adjacent_sides(v);
    return is_boundary(side1) || is_boundary(side2);
}

bool Block::is_boundary(Side side) const {
    return boundary(side) != nullptr;
}

void Block::set_boundary(Side side, Curve::Ref curve) {
    m_boundaries[static_cast<int>(side)] = curve;
    m_adjacent_blocks[static_cast<int>(side)].reset();
}

void Block::set_boundary(BaseNode::Ref v1, BaseNode::Ref v2, Curve::Ref curve) {
    set_boundary(get_side(v1, v2), curve);
}

// ---------- Выбор сторон и смежные блоки ----------- -----------------------------------------------------------------

Side Block::get_side(BaseNode::Ref v1, BaseNode::Ref v2) const {
    // Глупый полный перебор, хотя тут всего по 4 узла/грани
    int idx_v1 = base_node_index(v1);
    int idx_v2 = base_node_index(v2);
    if (idx_v1 == idx_v2) {
        throw std::runtime_error("Invalid side, same nodes");
    }
    if (idx_v2 < idx_v1) {
        std::swap(idx_v1, idx_v2);
    }

    for (Side side: sides_2D) {
        int idx_w1 = node_idx(side, 0);
        int idx_w2 = node_idx(side, 1);
        if (idx_w2 < idx_w1) {
            std::swap(idx_w1, idx_w2);
        }

        if (idx_v1 == idx_w1 && idx_v2 == idx_w2) {
            return side;
        }
    }

    for (Side side: sides_2D) {
        int idx_w1 = node_idx(side, 0);
        int idx_w2 = node_idx(side, 1);
        if (idx_w2 < idx_w1) {
            std::swap(idx_w1, idx_w2);
        }

        if (idx_v1 == idx_w1 && idx_v2 == idx_w2) {
            return side;
        }
    }
    throw std::runtime_error("Can't find face between two vertices");
}

Axis Block::get_axis(BaseNode::Ref v1, BaseNode::Ref v2) const {
    return to_axis(get_side(v1, v2));
}

Block::Ptr Block::adjacent_block(Side side) const {
    return m_adjacent_blocks[static_cast<int>(side)].lock();
}

Side Block::twin_face(Side side) const {
    int f_idx = side_to_idx(side);
    return idx_to_side((f_idx - m_rotations[static_cast<int>(side)] + 4) % 4);
}

// ---------- Внутренние вершины блока ---------------------------------------------------------------------------------

BsVertex::Ptr& Block::corner_vertex(int v_idx) {
    switch (v_idx) {
        case 0: return vertex( 0,  0);
        case 1: return vertex(-1,  0);
        case 2: return vertex( 0, -1);
        case 3: return vertex(-1, -1);
        default: throw std::runtime_error("Invalid base node index");
    }
}

BsVertex::Ptr& Block::boundary_vertex(Side side, int idx) {
    switch (side) {
        case Side::B: return vertex(idx, 0);
        case Side::R: return vertex(-1, idx);
        case Side::T: return vertex(-idx - 1, -1);
        case Side::L: return vertex(0, -idx - 1);
        default: throw std::runtime_error("Invalid base node index");
    }
}

BsVertex::Ptr& Block::near_boundary_vertex(Side side, int idx) {
    switch (side) {
        case Side::B: return vertex(idx, 1);
        case Side::R: return vertex(-2, idx);
        case Side::T: return vertex(-idx - 1, -2);
        case Side::L: return vertex(1, -idx - 1);
        default: throw std::runtime_error("Invalid base node index");
    }
}

// ---------- Конформные приколы ---------------------------------------------------------------------------------------

double Block::modulus() const { return m_modulus; }

void Block::estimate_modulus() {
    double len1 = length(Side::B) + length(Side::T);
    double len2 = length(Side::L) + length(Side::R);
    set_modulus(len1 / len2);
}

inline double restricted_modulus(double K) {
    static constexpr double max_modulus = 100.0;
    static constexpr double min_modulus = 1.0 / max_modulus;
    if (K > max_modulus) {
        std::cerr << "Modulus exceeded max\n";
        K = max_modulus;
    }
    if (K < min_modulus) {
        std::cerr << "Modulus exceeded min\n";
        K = min_modulus;
    }
    return K;
}

void Block::set_modulus(double K) {
    if (bad_number(K)) {
        throw std::runtime_error("Block::update_ratio: bad modulus");
    }
    m_modulus = restricted_modulus(K);
    if (size1() >= 1 && size2() >= 1) {
        m_lambda[1] = m_modulus * size2() / size1();
        m_lambda[0] = 1.0 / m_lambda[1];
    }
    else {
        m_lambda[0] = m_lambda[1] = 0.0;
    }
}

double* Block::lambda_ptr(Axis axis) { return &m_lambda[static_cast<int>(axis)]; }

double* Block::lambda_ptr(Side side) { return &m_lambda[static_cast<int>(to_axis(side))]; }

void Block::update_lambda() {
    if (good_number(m_modulus) && size1() >= 1 && size2() >= 1) {
        m_lambda[1] = m_modulus * size2() / size1();
        m_lambda[0] = 1.0 / m_lambda[1];
    }
    else {
        m_lambda[0] = m_lambda[1] = NAN;
    }
}

std::string Block::conformal_info() const {
    return std::format("Block {:3}.  K: {:.2f};  λ_1: {:.2f}; λ_2: {:.2f}", m_index, m_modulus, m_lambda[0], m_lambda[1]);
}

// ---------- Функции "верхнего уровня" --------------------------------------------------------------------------------

void Block::clear_vertices() {
    set_size(Axis::X, 0);
    set_size(Axis::Y, 0);
    m_vertices = {};
}

void Block::link(Block::Ref B1, Block::Ref B2) {
    if (!B1 || !B2) {
        throw std::runtime_error("Invalid block");
    }
    if (B1 == B2) {
        throw std::runtime_error("Can't link same blocks");
    }

    // Ищем просто полным перебором
    for (Side side1: sides_2D) {
        BaseNode::Ptr a1 = B1->base_node(node_idx(side1, 0));
        BaseNode::Ptr b1 = B1->base_node(node_idx(side1, 1));
        if (a1.get() > b1.get()) std::swap(a1, b1);

        for (Side side2: sides_2D) {
            BaseNode::Ptr a2 = B2->base_node(node_idx(side2, 0));
            BaseNode::Ptr b2 = B2->base_node(node_idx(side2, 1));
            if (a2.get() > b2.get()) std::swap(a2, b2);

            if (a1 == a2 && b1 == b2) {
                if (B1->boundary(side1)) {
                    throw std::runtime_error("Try link through boundary #1");
                }
                B1->m_adjacent_blocks[static_cast<int>(side1)] = B2;
                B1->m_rotations[static_cast<int>(side1)] = rotation(side1, side2);

                if (B2->boundary(side2)) {
                    throw std::runtime_error("Try link through boundary #2");
                }
                B2->m_adjacent_blocks[static_cast<int>(side2)] = B1;
                B2->m_rotations[static_cast<int>(side2)] = rotation(side2, side1);
            }
        }
    }
}

void Block::create_vertices_init() {
    using geom::SqQuad;

    check_consistency();
    if (size1() * size2() > 1000000) {
        throw std::runtime_error("Too much vertices");
    }

    m_vertices.resize({size1() + 1, size2() + 1}, nullptr);

    Vector3d v0 = base_node(0)->pos();
    Vector3d v1 = base_node(1)->pos();
    Vector3d v2 = base_node(2)->pos();
    Vector3d v3 = base_node(3)->pos();
    Vector3d vL = center(Side::L);
    Vector3d vR = center(Side::R);
    Vector3d vB = center(Side::B);
    Vector3d vT = center(Side::T);

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
        for (Side side: sides_2D) {
            if (!boundary(side)) { continue; }

            for (int idx = 1; idx < size(side); ++idx) {
                Vector3d v = (boundary_vertex(side, idx - 1)->pos +
                              boundary_vertex(side, idx + 1)->pos) / 2.0;
                v = boundary(side)->projection(v);
                boundary_vertex(side, idx)->next = v;
            }
        }

        for (Side side: sides_2D) {
            if (!boundary(side)) { continue; }
            for (int idx = 1; idx < size(side); ++idx) {
                boundary_vertex(side, idx)->pos = boundary_vertex(side, idx)->next;
            }
        }
    }

    // Генерация внутренних вершин, Coons patch.
    for (int i = 1; i < size1(); ++i) {
        for (int j = 1; j < size2(); ++j) {
            double t = double(i) / size1();
            double s = double(j) / size2();

            Vector3d Lc = (1 - t) * vertex(0, j)->pos + t * vertex(-1, j)->pos;
            Vector3d Ld = (1 - s) * vertex(i, 0)->pos + s * vertex(i, -1)->pos;
            Vector3d B = (1 - t) * (1 - s) * vertex(0, 0)->pos + (1 - t) * s * vertex(0, -1)->pos +
                         t * (1 - s) * vertex(-1, 0)->pos + t * s * vertex(-1, -1)->pos;
            Vector3d v = Lc + Ld - B;
            vertex(i, j) = BsVertex::create(v.x(), v.y());
        }
    }
}

void Block::create_vertices_again() {
    check_consistency();
    if (size1() * size2() > 1000000) {
        throw std::runtime_error("Too much vertices");
    }

    if (m_vertices.empty()) {
        throw std::runtime_error("Empty previous vertices array");
    }

    if (m_vertices.size1() < 2 || m_vertices.size2() < 2) {
        throw std::runtime_error("Empty previous vertices array");
    }

    // Берем прошлые вершины за основу
    Array2D<Vector3d> prev(m_vertices.sizes(), Vector3d::Zero());
    for (int i = 0; i < m_vertices.size1(); ++i) {
        for (int j = 0; j < m_vertices.size2(); ++j) {
            prev(i, j) = m_vertices(i, j)->pos;
        }
    }

    auto bilinear = [&prev](double x, double y) -> Vector3d {
        double i = x * (prev.size1() - 1);
        double j = y * (prev.size2() - 1);
        int i1 = std::min(static_cast<int>(std::floor(i)), prev.size1() - 2);
        int j1 = std::min(static_cast<int>(std::floor(j)), prev.size2() - 2);
        int i2 = i1 + 1;
        int j2 = j1 + 1;
        return (j2 - j) * ((i2 - i) * prev(i1, j1) + (i - i1) * prev(i2, j1)) +
               (j - j1) * ((i2 - i) * prev(i1, j2) + (i - i1) * prev(i2, j2));
    };

    m_vertices = Array2D<BsVertex::Ptr>({size1() + 1, size2() + 1}, nullptr);
    for (int i = 0; i <= size1(); ++i) {
        for (int j = 0; j <= size2(); ++j) {
            double x = math::between(double(i) / size1(), 0.0, 1.0);
            double y = math::between(double(j) / size2(), 0.0, 1.0);
            m_vertices(i,  j) = BsVertex::create(bilinear(x, y));
        }
    }

    // Проекция на границах
    for (Side side: sides_2D) {
        if (!boundary(side)) { continue; }

        for (int idx = 0; idx <= size(side); ++idx) {
            Vector3d v = boundary_vertex(side, idx)->pos;
            boundary_vertex(side, idx)->pos = boundary(side)->projection(v);
        }
    }
}

void Block::merge_vertices() {
    // Связать вершины с соседом
    for (Side side: sides_2D) {
        auto neib = adjacent_block(side);

        if (!neib) { continue; }

        // Сосед есть, но без вершин
        if (neib->vertices().empty()) {
            continue;
        }

        Side twin = twin_face(side);

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
            std::vector edges {
                BsEdge::Inside(vertex(i + 1, j), lambda_ptr(Axis::X)),
                BsEdge::Inside(vertex(i, j + 1), lambda_ptr(Axis::Y)),
                BsEdge::Inside(vertex(i - 1, j), lambda_ptr(Axis::X)),
                BsEdge::Inside(vertex(i, j - 1), lambda_ptr(Axis::Y))
            };
            vertex(i, j)->set_edges(edges);
        }
    }

    // Вершины на границах блока (без угловых)
    for (Side side: sides_2D) {
        if (!boundary(side) && !adjacent_block(side)) {
            throw std::runtime_error("BFace is not boundary and has no neighbor");
        }

        Axis axis = to_axis(side);
        Axis perp_axis = opposite(axis);

        if (boundary(side)) {
            for (int idx = 1; idx < size(side); ++idx) {
                // Добавить границу
                boundary_vertex(side, idx)->add_boundary(boundary(side).get());
                // Добавить связи (три штуки)
                std::vector edges {
                    BsEdge::Border(boundary_vertex(side, idx + 1), lambda_ptr(axis)),
                    BsEdge::Inside(near_boundary_vertex(side, idx),  lambda_ptr(perp_axis)),
                    BsEdge::Border(boundary_vertex(side, idx - 1), lambda_ptr(axis))
                };
                boundary_vertex(side, idx)->set_edges(edges);
            }
        } else {
            auto neib = adjacent_block(side);
            if (!neib) {
                throw std::runtime_error("Face is not boundary and has no neighbor");
            }

            int N = size(side);
            Side twin = twin_face(side);

            Axis neib_axis = to_axis(twin);
            Axis prep_neib_axis = opposite(neib_axis);
            for (int idx = 1; idx < N; ++idx) {
                std::vector edges {
                    BsEdge::Inside(boundary_vertex(side, idx + 1), lambda_ptr(axis), neib->lambda_ptr(neib_axis)),
                    BsEdge::Inside(near_boundary_vertex(side, idx), lambda_ptr(perp_axis)),
                    BsEdge::Inside(boundary_vertex(side, idx - 1), lambda_ptr(axis), neib->lambda_ptr(neib_axis)),
                    BsEdge::Inside(neib->near_boundary_vertex(twin, N - idx), neib->lambda_ptr(prep_neib_axis))
                };
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
            boundary(side1) && boundary(side2)) {
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
                    block->boundary_vertex(side, 1),
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
                    block->boundary_vertex(side, 1),
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
                    block->boundary_vertex(side, block->size(side) - 1),
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

            auto corner = corner_vertex(v_idx);
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
                    block->boundary_vertex(side, 1),
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

void Block::check_consistency() const {
    if (size1() < 1 || size2() < 1) {
        throw std::runtime_error("Zero size of some block");
    }
    for (Side side: sides_2D) {
        if (boundary(side)) { continue; }

        auto neib = adjacent_block(side);
        if (neib) {
            Side twin = twin_face(side);

            int size1 = size(side);
            int size2 = neib->size(twin);

            if (size1 != size2) {
                std::string message = std::format("Block::check_consistency error: size mismatch (blocks {}, {})", m_index, neib->index());
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

void Block::smooth() {
    for (int i = 0; i <= size1(); ++i) {
        for (int j = 0; j <= size2(); ++j) {
            auto vc = vertex(i, j);
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
                Curve* boundary = vc->boundary();
                Vector3d next = Vector3d::Zero();
                double sum = 0.0;
                for (auto edge: vertex(i, j)->adjacent()) {
                    if (edge.boundary()) {
                        next += edge.lambda() * edge.pos();
                        sum += edge.lambda();
                    }
                    else {
                        next += 2.0 * edge.lambda() * boundary->projection(edge.pos());
                        sum += 2.0 * edge.lambda();
                    }
                }
                next /= sum;
                vc->next = boundary->projection(next);
                continue;
            }

            // Внутренняя вершина
            Vector3d next = Vector3d::Zero();
            double sum = 0.0;
            for (auto edge: vertex(i, j)->adjacent()) {
                next += edge.lambda() * edge.pos();
                sum += edge.lambda();
            }
            next /= sum;
            vc->next = next;
        }
    }
}

double Block::move_vertices() {
    z_assert(m_vertices.size1() == size1() + 1, "Bad size1");
    z_assert(m_vertices.size2() == size2() + 1, "Bad size2");

    double err = 0.0;

    for (int i = 0; i <= size1(); ++i) {
        for (int j = 0; j <= size2(); ++j) {
            BsVertex::Ref vc = vertex(i, j);

            double L = 0.0;
            for (const auto edge: vertex(i, j)->adjacent()) {
                L = std::max(L, (vc->pos - edge.pos()).norm());
            }
            err = std::max(err, (vc->pos - vc->next).norm() / L);

            vc->pos = vc->next;
        }
    }
    return err;
}

void Block::update_modulus() {
    int Nx = size1();
    int Ny = size2();

    double sum_x = 0.0;
    double sum_y = 0.0;
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            Vector3d v1 = vertex(i, j)->pos;
            Vector3d v2 = vertex(i+1, j)->pos;
            Vector3d v3 = vertex(i+1, j+1)->pos;
            Vector3d v4 = vertex(i, j+1)->pos;

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

    double K = std::sqrt(sum_y / sum_x);
    set_modulus(K);
}

} // namespace zephyr::geom::generator
