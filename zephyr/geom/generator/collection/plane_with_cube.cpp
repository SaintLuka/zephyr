#include <iostream>
#include <algorithm>
#include <cassert>
#include <memory>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>

#include <zephyr/geom/generator/curve/plane.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/collection/plane_with_cube.h>

namespace zephyr::geom::generator::collection {

PlaneWithCube::PlaneWithCube(
        double xmin, double xmax,
        double ymin, double ymax,
        double xc, double yc,
        double r) :
        Generator("collection.plane-with-cube"),
        m_xmin(xmin), m_xmax(xmax),
        m_ymin(ymin), m_ymax(ymax),
        m_xc(xc), m_yc(yc),
        m_r(r) {

    double a = r;
    double b = 1.7 * a;
    
    double x1 = m_xmin;
    double x2 = m_xc - b;
    double x3 = m_xc;
    double x4 = m_xc + b;
    double x5 = m_xmax;
    double y1 = m_ymin;
    double y2 = m_yc - b;
    double y3 = m_yc;
    double y4 = m_yc + b;
    double y5 = m_ymax;

    double d = 0.5 * a;

    v01 = BaseNode::create(x1, y1);
    v02 = BaseNode::create(x2, y1);
    v03 = BaseNode::create(x3, y1);
    v04 = BaseNode::create(x4, y1);
    v05 = BaseNode::create(x5, y1);
    v06 = BaseNode::create(x1, y2);
    v07 = BaseNode::create(x2, y2);
    v08 = BaseNode::create(x3, y2);
    v09 = BaseNode::create(x4, y2);
    v10 = BaseNode::create(x5, y2);
    v11 = BaseNode::create(x1, y3);
    v12 = BaseNode::create(x2, y3);
    v13 = BaseNode::create(x3-a, y3);
    v14 = BaseNode::create(x3-d, y3-d);
    v15 = BaseNode::create(x3, y3-a);
    v16 = BaseNode::create(x3+d, y3-d);
    v17 = BaseNode::create(x3+a, y3);
    v18 = BaseNode::create(x4, y3);
    v19 = BaseNode::create(x5, y3);
    v20 = BaseNode::create(x1, y4);
    v21 = BaseNode::create(x2, y4);
    v22 = BaseNode::create(x3-d, y3+d);
    v23 = BaseNode::create(x3, y3+a);
    v24 = BaseNode::create(x3, y4);
    v25 = BaseNode::create(x3+d, y3+d);
    v26 = BaseNode::create(x4, y4);
    v27 = BaseNode::create(x5, y4);
    v28 = BaseNode::create(x1, y5);
    v29 = BaseNode::create(x2, y5);
    v30 = BaseNode::create(x3, y5);
    v31 = BaseNode::create(x4, y5);
    v32 = BaseNode::create(x5, y5);

    m_blocks += {v01, v02, v06, v07};
    m_blocks += {v02, v03, v07, v08};
    m_blocks += {v03, v04, v08, v09};
    m_blocks += {v04, v05, v09, v10};

    m_blocks += {v06, v07, v11, v12};
    m_blocks += {v07, v14, v12, v13};
    m_blocks += {v07, v08, v14, v15};
    m_blocks += {v08, v09, v15, v16};
    m_blocks += {v16, v09, v17, v18};
    m_blocks += {v09, v10, v18, v19};

    m_blocks += {v11, v12, v20, v21};
    m_blocks += {v12, v13, v21, v22};
    m_blocks += {v22, v21, v23, v24};
    m_blocks += {v23, v25, v24, v26};
    m_blocks += {v17, v18, v25, v26};
    m_blocks += {v18, v19, v26, v27};

    m_blocks += {v20, v21, v28, v29};
    m_blocks += {v21, v24, v29, v30};
    m_blocks += {v24, v26, v30, v31};
    m_blocks += {v26, v27, v31, v32};

    left  = Plane::create(v01, v28);
    right = Plane::create(v05, v32);
    bottom = Plane::create(v01, v05);
    top   = Plane::create(v28, v32);
    side1 = Plane::create(v13, v15);
    side2 = Plane::create(v15, v17);
    side3 = Plane::create(v13, v23);
    side4 = Plane::create(v23, v17);

    m_blocks.set_boundary({v01, v06, v11, v20, v28}, left);
    m_blocks.set_boundary({v05, v10, v19, v27, v32}, right);
    m_blocks.set_boundary({v01, v02, v03, v04, v05}, bottom);
    m_blocks.set_boundary({v28, v29, v30, v31, v32}, top);
    m_blocks.set_boundary({v13, v14, v15}, side1);
    m_blocks.set_boundary({v15, v16, v17}, side2);
    m_blocks.set_boundary({v13, v22, v23}, side3);
    m_blocks.set_boundary({v23, v25, v17}, side4);

    //m_blocks.plot_layout();
    //m_blocks.set_verbosity(5);
    m_blocks.set_iters_count(100);
    m_blocks.optimize();
}

void PlaneWithCube::set_nx(int Nx) {
    m_blocks.set_size({v01, v02, v03, v04, v05}, Nx); // Нижняя граница
    m_blocks.set_size({v28, v29, v30, v31, v32}, Nx); // Верхняя граница
}

void PlaneWithCube::set_ny(int Ny) {
    m_blocks.set_size({v01, v06, v11, v20, v28}, Ny); // Левая граница
    m_blocks.set_size({v05, v10, v19, v27, v32}, Ny); // Правая граница
}

void PlaneWithCube::set_boundaries(Boundaries bounds) const {
    left->set_boundary(bounds.left);
    right->set_boundary(bounds.right);
    bottom->set_boundary(bounds.bottom);
    top->set_boundary(bounds.top);
    side1->set_boundary(bounds.hole);
    side2->set_boundary(bounds.hole);
    side3->set_boundary(bounds.hole);
    side4->set_boundary(bounds.hole);
}

void PlaneWithCube::set_axial(bool axial) {
    m_blocks.set_axial(axial);
}

void PlaneWithCube::set_adaptive(bool adaptive) {
    m_blocks.set_adaptive(adaptive);
}

void PlaneWithCube::set_linear(bool linear) {
    m_blocks.set_linear(linear);
}

Box PlaneWithCube::bbox() const {
    Vector3d vmin(m_xmin, m_ymin, 0.0);
    Vector3d vmax(m_xmax, m_ymax, 0.0);
    return {vmin, vmax};
}

void PlaneWithCube::check_params() const {
    if (m_xmin >= m_xmax) {
        throw std::runtime_error("PlaneWithCube Error: x_min >= x_max");
    }
    if (m_ymin >= m_ymax) {
        throw std::runtime_error("PlaneWithCube Error: y_min >= y_max");
    }

    double R = 1.5 * m_r;
    if (m_xmin >= m_xc - R || m_xmax <= m_xc + R || m_ymin >= m_yc - R || m_ymax <= m_yc + R) {
        throw std::runtime_error("PlaneWithCube Error: big hole");
    }
}

Grid PlaneWithCube::make() const {
    return m_blocks.make();
}

} // namespace zephyr::geom::generator::collection