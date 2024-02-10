#include <algorithm>

#include <zephyr/geom/generator/curve/curve.h>
#include <zephyr/geom/generator/bs_vertex.h>


namespace zephyr::geom::generator {

BsVertex::BsVertex(const Vector3d &v)
    : index(-1), v1(v), v2(v) {
}

BsVertex::BsVertex(double x, double y)
    : index(-1), v1(x, y, 0.0), v2(x, y, 0.0) {
}

BsVertex::Ptr BsVertex::create(const Vector3d &v) {
    return std::make_shared<BsVertex>(v);
}

BsVertex::Ptr BsVertex::create(double x, double y) {
    return std::make_shared<BsVertex>(x, y);
}

int BsVertex::n_adjacent() const {
    return m_adjacent.size();
}

void BsVertex::fix() {
    m_adjacent.clear();
}

const std::vector<BsVertex *> &BsVertex::adjacent_vertices() const {
    return m_adjacent;
}

void BsVertex::set_adjacent_vertices(const std::vector<BsVertex::Ptr> &vertices) {
    std::vector<std::pair<double, BsVertex::Ptr>> vec_verts;
    vec_verts.reserve(vertices.size());
    for (auto &v: vertices) {
        vec_verts.emplace_back(std::make_pair(std::atan2(v->v1.y() - v1.y(), v->v1.x() - v1.x()), v));
    }

    for (size_t i = 1; i < vec_verts.size(); ++i) {
        if (vec_verts[i].first < vec_verts[0].first) {
            vec_verts[i].first += 2.0 * M_PI;
        }
    }

    std::sort(vec_verts.begin(), vec_verts.end(),
              [](const std::pair<double, BsVertex::Ptr> &p1, const
              std::pair<double, BsVertex::Ptr> &p2) {
                  return p1.first < p2.first;
              });

    m_adjacent.reserve(vec_verts.size());
    for (auto &p: vec_verts) {
        m_adjacent.emplace_back(p.second.get());
    }
}

bool BsVertex::inner() const {
    return m_boundaries.empty();
}

bool BsVertex::corner() const {
    return m_boundaries.size() > 1;
}

Curve *BsVertex::boundary() const {
    if (m_boundaries.size() == 1) {
        return *m_boundaries.begin();
    }
    else {
        throw std::runtime_error("BsVertex error: #1464");
    }
}

std::set<Boundary> BsVertex::boundaries() const {
    std::set<Boundary> res;
    for (auto& curve: m_boundaries) {
        res.insert(curve->boundary());
    }
    return res;
}

void BsVertex::add_boundary(Curve *boundary) {
    m_boundaries.insert(boundary);
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

} // namespace zephyr::geom::generator
