#include <random>

#include <zephyr/geom/vector.h>
#include <zephyr/mesh/decomp/rwalk.h>

using zephyr::geom::Vector3d;

namespace zephyr::mesh::decomp {

RWalk::RWalk(const Box &domain, int size)
    : Decomposition(size), m_domain(domain) {
    const int multiplier = 10;

    m_diagram = VDiagram(domain, multiplier * size);

    if (m_domain.is_2D()) {
        m_step = 0.3 * std::sqrt(m_domain.area() / m_diagram.size());
    } else {
        m_step = 0.3 * std::cbrt(m_domain.volume() / m_diagram.size());
    }
}

int RWalk::rank(const EuCell &elem) const {
    return m_diagram.rank(elem.center()) % m_size;
}

void RWalk::balancing(const std::vector<double> &w) {
    static std::default_random_engine gen;

    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    for (int i = 0; i < m_diagram.size(); ++i) {
        double x = m_diagram.get_coord_x(i);
        double y = m_diagram.get_coord_y(i);

        x += 2.0 * m_step * (uniform(gen) - 0.5);
        y += 2.0 * m_step * (uniform(gen) - 0.5);

        if (x < m_domain.vmin.x() || x > m_domain.vmax.x()) {
            x = m_domain.vmin.x() + m_domain.size().x() * uniform(gen);
        }
        if (y < m_domain.vmin.y() || y > m_domain.vmax.y()) {
            y = m_domain.vmin.y() + m_domain.size().y() * uniform(gen);
        }

        m_diagram.set_coords(i, x, y);
    }
}

} // namespace zephyr::mesh::decomp