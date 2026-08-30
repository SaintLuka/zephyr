#include <random>

#include <zephyr/geom/vector.h>
#include <zephyr/mesh/decomp/rwalk.h>

using zephyr::geom::Vector3d;

namespace zephyr::mesh::decomp {

RWalk::RWalk(const Box &domain, int size)
    : Decomposition(size), domain_(domain) {
    constexpr int multiplier = 10;

    diagram_ = VDiagram(domain, multiplier * size);

    if (domain_.is_2D()) {
        step_ = 0.3 * std::sqrt(domain_.area() / diagram_.size());
    } else {
        step_ = 0.3 * std::cbrt(domain_.volume() / diagram_.size());
    }
}

int RWalk::rank(const EuCell &elem) const {
    return diagram_.rank(elem.center()) % m_size;
}

void RWalk::balancing(const std::vector<double> &w) {
    static std::default_random_engine gen;

    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    for (int i = 0; i < diagram_.size(); ++i) {
        Vector3d p = diagram_.get_coord(i);

        p.x() += 2.0 * step_ * (uniform(gen) - 0.5);
        p.y() += 2.0 * step_ * (uniform(gen) - 0.5);
        p.z() += 2.0 * step_ * (uniform(gen) - 0.5);

        diagram_.set_coords(i, p);
    }
}

} // namespace zephyr::mesh::decomp