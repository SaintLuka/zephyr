#include <iostream>
#include <algorithm>
#include <cassert>
#include <memory>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>

#include <zephyr/geom/generator/curve/plane.h>
#include <zephyr/geom/generator/collection/wedge.h>

namespace zephyr::geom::generator::collection {

Wedge::Wedge(
        double xmin, double xmax,
        double ymin, double ymax,
        double xw, double phi) :
        Generator("collection.wedge"),
        m_xmin(xmin), m_xmax(xmax),
        m_ymin(ymin), m_ymax(ymax),
        m_xw(xw), m_phi(phi) {

    check_params();

    double m_yw = m_ymin + (m_xmax - m_xw) * std::tan(m_phi);

    // Задаем базисные вершины для структурированных блоков
    v1 = BaseNode::create(m_xmin, m_ymin);
    v2 = BaseNode::create(m_xw,   m_ymin);
    v3 = BaseNode::create(m_xmax, m_yw);
    v4 = BaseNode::create(m_xmin, m_ymax);
    v5 = BaseNode::create(m_xw,   m_ymax);
    v6 = BaseNode::create(m_xmax, m_ymax);

    // Ограничивающие прямые области
    left   = Plane::create(v1, v4);
    right  = Plane::create(v3, v6);
    bottom = Plane::create(v1, v2);
    top    = Plane::create(v4, v6);
    wedge  = Plane::create(v2, v3);

    // Генератор сетки
    m_blocks += {v1, v2, v4, v5};
    m_blocks += {v2, v3, v5, v6};

    m_blocks.set_boundary(v1, v4, left);
    m_blocks.set_boundary(v3, v6, right);
    m_blocks.set_boundary(v1, v2, bottom);
    m_blocks.set_boundary(v2, v3, wedge);
    m_blocks.set_boundary({v4, v5, v6}, top);

    m_blocks.set_verbosity(0);
    m_blocks.optimize({.steps=2});
}

void Wedge::check_params() const {
    if (m_xmin >= m_xmax) {
        throw std::runtime_error("Wedge::check_params: x_min >= x_max");
    }
    if (m_ymin >= m_ymax) {
        throw std::runtime_error("Wedge::check_params: y_min >= y_max");
    }
    if (m_xw >= m_xmax || m_xw <= m_xmin) {
        throw std::runtime_error("Wedge::check_params: x_w not in [x_min, x_max]");
    }
    if (m_phi < 0.0 || m_phi > 0.75 * M_PI) {
        throw std::runtime_error("Wedge::check_params: phi not in [0.0, pi/2]");
    }

    // высота клина
    double H = (m_xmax - m_xw) * std::tan(m_phi);
    if (H > 0.9 * (m_ymax - m_ymin)) {
        throw std::runtime_error("Wedge::check_params: big wedge, increase (y_max - y_min) or decrease angle");
    }
}

void Wedge::set_nx(int Nx) {
    m_blocks.set_size({v4, v5, v6}, Nx);
}

void Wedge::set_ny(int Ny) {
    m_blocks.set_size(v1, v4, Ny);
}

void Wedge::set_boundaries(Boundaries bounds) const {
    left->set_boundary(bounds.left);
    right->set_boundary(bounds.right);
    bottom->set_boundary(bounds.bottom);
    top->set_boundary(bounds.top);
    wedge->set_boundary(bounds.bottom);
}

Box Wedge::bbox() const {   
    Vector3d vmin(m_xmin, m_ymin, 0.0);
    Vector3d vmax(m_xmax, m_ymax, 0.0);

    return {vmin, vmax};
}

} // namespace zephyr::geom::generator::collection