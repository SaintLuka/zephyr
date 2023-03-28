#include <algorithm>

#include <zephyr/mesh/generator/vertex.h>


namespace zephyr { namespace mesh { namespace generator {

Vertex::Vertex(const Vector3d &v)
    : v1(v), v2(v), m_boundary(nullptr) {
}

Vertex::Vertex(double x, double y)
    : v1(x, y, 0.0), v2(x, y, 0.0), m_boundary(nullptr) {
}

Vertex::Ptr Vertex::create(const Vector3d &v) {
    return std::make_shared<Vertex>(v);
}

Vertex::Ptr Vertex::create(double x, double y) {
    return std::make_shared<Vertex>(x, y);
}

int Vertex::n_adjacent() const {
    return m_adjacent.size();
}

const std::vector<Vertex *> &Vertex::adjacent_vertices() const {
    return m_adjacent;
}

void Vertex::set_adjacent_vertices(const std::vector<Vertex::Ptr> &vertices) {
    std::vector<std::pair<double, Vertex::Ptr>> vec_verts;
    vec_verts.reserve(vertices.size());
    for (auto &v: vertices) {
        vec_verts.emplace_back(std::make_pair(std::atan2(v->v1.y() - v1.y(), v->v1.x() - v1.x()), v));
    }

    for (int i = 1; i < vec_verts.size(); ++i) {
        if (vec_verts[i].first < vec_verts[0].first) {
            vec_verts[i].first += 2.0 * M_PI;
        }
    }

    std::sort(vec_verts.begin(), vec_verts.end(),
              [](const std::pair<double, Vertex::Ptr> &p1, const
              std::pair<double, Vertex::Ptr> &p2) {
                  return p1.first < p2.first;
              });

    m_adjacent.reserve(vec_verts.size());
    for (auto &p: vec_verts) {
        m_adjacent.emplace_back(p.second.get());
    }
}

Curve *Vertex::boundary() const {
    return m_boundary;
}

void Vertex::set_boundary(Curve *boundary) {
    m_boundary = boundary;
}

BaseVertex::BaseVertex(const Vector3d &v, bool fixed)
    : m_v(v), m_fixed(fixed) {
}

BaseVertex::Ptr BaseVertex::create(const Vector3d &v, bool fixed) {
    return std::make_shared<BaseVertex>(v, fixed);
}

BaseVertex::Ptr BaseVertex::create(double x, double y, bool fixed) {
    return std::make_shared<BaseVertex>(Vector3d(x, y, 0.0), fixed);
}

bool BaseVertex::is_fixed() const {
    return m_fixed;
}

const Vector3d &BaseVertex::v() const {
    return m_v;
}

int BaseVertex::degree() const {
    return m_adjacent_blocks.size();
}

void BaseVertex::add_adjacent_block(Block *block) {
    for (auto b: m_adjacent_blocks) {
        if (b == block) {
            return;
        }
    }
    m_adjacent_blocks.push_back(block);
}

const std::vector<Block *> &BaseVertex::adjacent_blocks() const {
    return m_adjacent_blocks;
}

} // namespace generator
} // namespace mesh
} // namespace zephyr
