#include <algorithm>
#include <ranges>

#include <zephyr/geom/generator/bs_vertex.h>
#include <zephyr/geom/generator/curve/curve.h>

namespace zephyr::geom::generator {

BsEdge::BsEdge(BsVertex* v, double* lambda_1, double* lambda_2)
    : neib(v), lambda1(lambda_1), lambda2(lambda_2) { }

BsEdge BsEdge::Border(BsVertex_Ref v, double* lambda) {
    return BsEdge{v.get(), lambda, nullptr};
}

BsEdge BsEdge::Inside(BsVertex_Ref v, double* lambda) {
    return BsEdge{v.get(), lambda, lambda};
}

BsEdge BsEdge::Inside(BsVertex::Ref v, double* lambda1, double* lambda2) {
    return BsEdge{v.get(), lambda1, lambda2};
}

double BsEdge::lambda() const {
    z_assert(lambda1 != nullptr, "BsEdge::lambda: lambda1 is null");
    if (boundary()) { return (*lambda1); }
    return std::sqrt((*lambda1) * (*lambda2));
}

BsVertex::BsVertex(const Vector3d &v)
    : index(-1), pos(v), next(v) { }

BsVertex::BsVertex(double x, double y)
    : index(-1), pos(x, y, 0.0), next(x, y, 0.0) { }

BsVertex::Ptr BsVertex::create(const Vector3d &v) {
    return std::make_shared<BsVertex>(v);
}

BsVertex::Ptr BsVertex::create(double x, double y) {
    return std::make_shared<BsVertex>(x, y);
}

void BsVertex::set_edges(const std::vector<BsEdge> &edges) {
    std::vector<std::pair<double, BsEdge>> sorted;
    sorted.reserve(edges.size());
    for (auto &edge: edges) {
        sorted.emplace_back(std::atan2(edge.y() - pos.y(), edge.x() - pos.x()), edge);
    }

    for (size_t i = 1; i < sorted.size(); ++i) {
        if (sorted[i].first < sorted[0].first) {
            sorted[i].first += 2.0 * M_PI;
        }
    }

    std::ranges::sort(sorted,
                      [](const std::pair<double, BsEdge> &p1,
                         const std::pair<double, BsEdge> &p2) {
                          return p1.first < p2.first;
                      });

    m_edges.clear();
    m_edges.reserve(sorted.size());
    for (auto& val: sorted | std::views::values) {
        m_edges.emplace_back(val);
    }
}

void BsVertex::add_boundary(Curve* boundary) {
    if (boundary) {
        m_boundaries.insert(boundary);
    }
}

Curve* BsVertex::boundary() const {
    if (m_boundaries.size() == 1) {
        return *m_boundaries.begin();
    }
    throw std::runtime_error("BsVertex::boundary: more than one boundary");
}

std::set<Boundary> BsVertex::boundaries() const {
    std::set<Boundary> res;
    for (auto& curve: m_boundaries) {
        res.insert(curve->boundary());
    }
    return res;
}

} // namespace zephyr::geom::generator
