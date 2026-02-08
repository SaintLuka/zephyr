#include <algorithm>
#include <ranges>

#include <zephyr/geom/generator/curve/curve.h>
#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/block.h>


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
    std::vector<std::pair<double, BsVertex::Ptr>> sorted;
    sorted.reserve(vertices.size());
    for (auto &v: vertices) {
        sorted.emplace_back(std::atan2(v->v1.y() - v1.y(), v->v1.x() - v1.x()), v);
    }

    for (size_t i = 1; i < sorted.size(); ++i) {
        if (sorted[i].first < sorted[0].first) {
            sorted[i].first += 2.0 * M_PI;
        }
    }

    std::ranges::sort(sorted,
                      [](const std::pair<double, BsVertex::Ptr> &p1,
                         const std::pair<double, BsVertex::Ptr> &p2) {
                          return p1.first < p2.first;
                      });

    m_adjacent.clear();
    m_adjacent.reserve(sorted.size());
    for (auto& val : sorted | std::views::values) {
        m_adjacent.emplace_back(val.get());
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
    if (boundary) {
        m_boundaries.insert(boundary);
    }
}

} // namespace zephyr::geom::generator
