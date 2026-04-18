#include <iostream>
#include <algorithm>
#include <cassert>
#include <memory>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>

#include <zephyr/geom/generator/curve/plane.h>
#include <zephyr/geom/generator/curve/circle.h>
#include <zephyr/geom/generator/collection/semicircle_cutout.h>

namespace zephyr::geom::generator::collection {

SemicircleCutout::SemicircleCutout(
        double xmin, double xmax,
        double ymin, double ymax,
        double xc, double r) :
        Generator("collection.semicircle-cutout"),
        m_xmin(xmin), m_xmax(xmax),
        m_ymin(ymin), m_ymax(ymax),
        m_xc(xc), m_r(r) {
    check_params();

    double r_max = std::min({m_xmax - m_xc, m_xc - m_xmin, m_ymax - m_ymin});
    double R = std::min(1.2 * m_r, r_max);
    double phi1 = 0.75 * M_PI;
    double phi2 = 0.25 * M_PI;

    //double R = 1.1 * m_r;

    double x1 = m_xc - R;// * std::cos(phi1);
    double x2 = m_xc + m_r * std::cos(phi1);
    double x3 = m_xc + m_r * std::cos(phi2);
    double x4 = m_xc + R;// * std::cos(phi2);
    double y1 = m_ymin + m_r * std::sin(phi1);
    double y2 = m_ymin + R; //R * std::sin(phi1);

    // Задаем базисные вершины для структурированных блоков
    v1 = BaseNode::create(m_xmin,   m_ymin);
    v2 = BaseNode::create(m_xc-R,   m_ymin);
    v3 = BaseNode::create(m_xc-m_r, m_ymin);
    v4 = BaseNode::create(m_xc+m_r, m_ymin);
    v5 = BaseNode::create(m_xc+R,   m_ymin);
    v6 = BaseNode::create(m_xmax,   m_ymin);

    v7  = BaseNode::create(m_xmin, y2);
    v8  = BaseNode::create(x1, y2);
    v9  = BaseNode::create(x2, y1);
    v10 = BaseNode::create(x3, y1);
    v11 = BaseNode::create(x4, y2);
    v12 = BaseNode::create(m_xmax, y2);

    v13 = BaseNode::create(m_xmin, m_ymax);
    v14 = BaseNode::create(x1, m_ymax);
    v15 = BaseNode::create(x4, m_ymax);
    v16 = BaseNode::create(m_xmax, m_ymax);

    // Кривые на границе
    left   = Plane::create(v1, v13);
    right  = Plane::create(v6, v16);
    bottom = Plane::create(v1, v6);
    top    = Plane::create(v13, v16);
    circle = Circle::create(v3, v9, v4);

    // Генератор сетки
    m_blocks += {v1, v2, v7, v8};
    m_blocks += {v2, v3, v8, v9};
    m_blocks += {v9, v10, v8, v11};
    m_blocks += {v4, v5, v10, v11};
    m_blocks += {v5, v6, v11, v12};
    m_blocks += {v7, v8, v13, v14};
    m_blocks += {v8, v11, v14, v15};
    m_blocks += {v11, v12, v15, v16};

    m_blocks.set_boundary({v1, v7, v13}, left);
    m_blocks.set_boundary({v6, v12, v16}, right);
    m_blocks.set_boundary({v1, v2, v3}, bottom);
    m_blocks.set_boundary({v4, v5, v6}, bottom);
    m_blocks.set_boundary({v3, v9, v10, v4}, circle);
    m_blocks.set_boundary({v13, v14, v15, v16}, top);

    m_blocks.set_verbosity(0);
    m_blocks.optimize({.steps=2});
}

void SemicircleCutout::check_params() const {
    if (m_xmin >= m_xmax) {
        throw std::runtime_error("SemicircleCutout::check_params: x_min >= x_max");
    }
    if (m_ymin >= m_ymax) {
        throw std::runtime_error("SemicircleCutout::check_params: y_min >= y_max");
    }
    if (m_xc >= m_xmax || m_xc <= m_xmin) {
        throw std::runtime_error("SemicircleCutout::check_params: x_c not in [x_min, x_max]");
    }
    if (1.5 * m_r > (m_ymax - m_ymin)) {
        // По оси y умещается 1.5 радиуса
        throw std::runtime_error("SemicircleCutout::check_params: big radius, increase (y_max - y_min)");
    }
    if (1.5 * m_r > (m_xmax - m_xc)) {
        // По оси x умещается 1.5 радиуса вправо
        throw std::runtime_error("SemicircleCutout::check_params: big radius, increase (x_max - x_c)");
    }
    if (1.5 * m_r > (m_xc - m_xmin)) {
        // По оси x умещается 1.5 радиуса влево
        throw std::runtime_error("SemicircleCutout::check_params: big radius, increase (x_c - x_min)");
    }
}

void SemicircleCutout::set_nx(int Nx) {
    m_blocks.set_size({v13, v14, v15, v16}, Nx);
}

void SemicircleCutout::set_ny(int Ny) {
    m_blocks.set_size({v1, v7, v13}, Ny);
}

void SemicircleCutout::set_boundaries(Boundaries bounds) const {
    left->set_boundary(bounds.left);
    right->set_boundary(bounds.right);
    bottom->set_boundary(bounds.bottom);
    top->set_boundary(bounds.top);
    circle->set_boundary(bounds.bottom);
}

void SemicircleCutout::set_axial(bool axial) {
    m_blocks.set_axial(axial);
}

void SemicircleCutout::set_adaptive(bool adaptive) {
    m_blocks.set_adaptive(adaptive);
}

void SemicircleCutout::set_linear(bool linear) {
    m_blocks.set_linear(linear);
}

Box SemicircleCutout::bbox() const {
    Vector3d vmin(m_xmin, m_ymin, 0.0);
    Vector3d vmax(m_xmax, m_ymax, 0.0);
    return {vmin, vmax};
}

Grid SemicircleCutout::make() const {
    return m_blocks.make();
}

} // namespace zephyr::geom::generator::collection