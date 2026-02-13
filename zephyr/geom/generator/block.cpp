#include <algorithm>
#include <array>
#include <memory>
#include <set>
#include <format>

#include <zephyr/geom/generator/block.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>

#include <zephyr/geom/intersection.h>
#include <zephyr/geom/primitives/triangle.h>
#include <zephyr/math/funcs.h>

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
        case Side::LEFT:   return std::array{2, 0}[idx];
        case Side::RIGHT:  return std::array{1, 3}[idx];
        case Side::BOTTOM: return std::array{0, 1}[idx];
        case Side::TOP:    return std::array{3, 2}[idx];
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
    return m_sizes[axis];
}

int Block::size(Side side) const {
    return m_sizes[side];
}

int Block::size(BaseNode::Ref v1, BaseNode::Ref v2) const {
    return m_sizes[get_side(v1, v2)];
}

void Block::set_size(Axis axis, int N) {
    if (m_sizes[axis] == N) return;
    m_sizes[axis] = N;
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
    return m_rel_sizes[axis];
}

double Block::rel_size(BaseNode::Ref v1, BaseNode::Ref v2) const {
    return m_rel_sizes[get_side(v1, v2)];
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
    if (m_rel_sizes[axis] == N) return;
    m_rel_sizes[axis] = N;
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

double* Block::lambda_ptr(Axis axis) { return &m_lambda[axis]; }

double* Block::lambda_ptr(Side side) { return &m_lambda[side]; }

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

Array2D<BsVertex::Ptr> Block::create_vertices(axis_pair<int> sizes) const {
    if (sizes[0] * sizes[1] > 1000000) {
        throw std::runtime_error("Too much vertices");
    }

    Array2D<BsVertex::Ptr> vertices({sizes[Axis::X] + 1, sizes[Axis::Y] + 1}, nullptr);

    for (Side side: sides_2D) {
        Vector3d p1 = base_node(node_idx(side, 0))->pos();
        Vector3d p2 = base_node(node_idx(side, 1))->pos();

        for (int i = 0; i <= sizes[side]; ++i) {
            double x = static_cast<double>(i) / sizes[side];
            vertices.boundary(side, i) = BsVertex::create((1.0 - x) * p1 + x * p2);
        }

        if (is_boundary(side)) {
            auto curve = boundary(side);
            for (int idx = 1; idx < sizes[side]; ++idx) {
                Vector3d p = vertices.boundary(side, idx)->pos;
                vertices.boundary(side, idx)->pos = curve->projection(p);
            }

            // Сглаживание границы
            for (int k = 0; k < 10; ++k) {
                for (int idx = 1; idx < sizes[side]; ++idx) {
                    Vector3d v = 0.5 * (vertices.boundary(side, idx - 1)->pos +
                                        vertices.boundary(side, idx + 1)->pos);
                    vertices.boundary(side, idx)->next = curve->projection(v);
                }

                for (int idx = 1; idx < sizes[side]; ++idx) {
                    vertices.boundary(side, idx)->pos = vertices.boundary(side, idx)->next;
                }
            }
        }
    }

    // Генерация внутренних вершин, Coons patch.
    for (int i = 1; i < sizes[Axis::X]; ++i) {
        for (int j = 1; j < sizes[Axis::Y]; ++j) {
            double t = static_cast<double>(i) / sizes[Axis::X];
            double s = static_cast<double>(j) / sizes[Axis::Y];

            Vector3d Lc = (1 - t) * vertices(0, j)->pos + t * vertices(-1, j)->pos;
            Vector3d Ld = (1 - s) * vertices(i, 0)->pos + s * vertices(i, -1)->pos;
            Vector3d B = (1 - t) * (1 - s) * vertices(0, 0)->pos + (1 - t) * s * vertices(0, -1)->pos +
                         t * (1 - s) * vertices(-1, 0)->pos + t * s * vertices(-1, -1)->pos;
            Vector3d v = Lc + Ld - B;
            vertices(i, j) = BsVertex::create(v.x(), v.y());
        }
    }
    return vertices;
}

Array2D<BsVertex::Ptr> Block::create_vertices_again(axis_pair<int> sizes) const {
    if (sizes[Axis::X] * sizes[Axis::Y] > 1000000) {
        throw std::runtime_error("Too much vertices");
    }

    if (m_mapping.empty()) {
        throw std::runtime_error("Empty previous vertices array");
    }

    if (m_mapping.size1() < 2 || m_mapping.size2() < 2) {
        throw std::runtime_error("Empty previous vertices array");
    }

    // Берем прошлые вершины за основу
    auto bilinear = [&prev=m_mapping](double x, double y) -> Vector3d {
        double i = x * (prev.size1() - 1);
        double j = y * (prev.size2() - 1);
        int i1 = std::min(static_cast<int>(std::floor(i)), prev.size1() - 2);
        int j1 = std::min(static_cast<int>(std::floor(j)), prev.size2() - 2);
        int i2 = i1 + 1;
        int j2 = j1 + 1;
        return (j2 - j) * ((i2 - i) * prev(i1, j1) + (i - i1) * prev(i2, j1)) +
               (j - j1) * ((i2 - i) * prev(i1, j2) + (i - i1) * prev(i2, j2));
    };

    Array2D<BsVertex::Ptr> vertices({sizes[Axis::X] + 1, sizes[Axis::Y] + 1}, nullptr);
    for (int i = 0; i <= sizes[Axis::X]; ++i) {
        for (int j = 0; j <= sizes[Axis::Y]; ++j) {
            double x = math::between(static_cast<double>(i) / sizes[Axis::X], 0.0, 1.0);
            double y = math::between(static_cast<double>(j) / sizes[Axis::Y], 0.0, 1.0);
            vertices(i,  j) = BsVertex::create(bilinear(x, y));
        }
    }

    // Проекция на границах
    for (Side side: sides_2D) {
        if (!boundary(side)) { continue; }

        auto curve = boundary(side);
        for (int idx = 0; idx <= sizes[side]; ++idx) {
            Vector3d v = vertices.boundary(side, idx)->pos;
            vertices.boundary(side, idx)->pos = curve->projection(v);
        }
    }
    return vertices;
}

void Block::update_modulus(const Array2D<BsVertex::Ptr>& vertices) {
    const int Nx = vertices.size(Axis::X) - 1;
    const int Ny = vertices.size(Axis::Y) - 1;

    double sum_x = 0.0;
    double sum_y = 0.0;
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            Vector3d v1 = vertices(i, j)->pos;
            Vector3d v2 = vertices(i+1, j)->pos;
            Vector3d v3 = vertices(i+1, j+1)->pos;
            Vector3d v4 = vertices(i, j+1)->pos;

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

void Block::set_mapping(const Array2D<BsVertex::Ptr>& vertices) {
    m_mapping.resize(vertices.sizes());
    for (int i = 0; i < vertices.size1(); ++i) {
        for (int j = 0; j < vertices.size2(); ++j) {
            m_mapping(i, j) = vertices(i, j)->pos;
        }
    }
}

} // namespace zephyr::geom::generator
