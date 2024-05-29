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
        BlockStructured(10),
        m_xmin(xmin), m_xmax(xmax),
        m_ymin(ymin), m_ymax(ymax),
        m_xc(xc), m_yc(yc),
        m_r(r), m_xi(2.0) {

    m_name = "collection.plane-with-cube";

    init_blocks();
}

void PlaneWithCube::set_nx(int Nx) {

    // TODO
    m_blocks[0].set_size(v1, v8, 5);
    m_blocks[0].set_size(v1, v6, 5);

    m_blocks[1].set_size(v1, v2, 5);

    m_blocks[2].set_size(v2, v3, 5);
    m_blocks[2].set_size(v3, v7, 5);

    m_blocks[3].set_size(v3, v11, 5);

    m_blocks[4].set_size(v3, v4, 5);

    m_blocks[5].set_size(v11, v18, 5);

    m_blocks[6].set_size(v14, v18, 5);

    m_blocks[7].set_size(v17, v18, 5);

    m_blocks[8].set_size(v16, v17, 5);
    m_blocks[8].set_size(v16, v13, 5);

    m_blocks[9].set_size(v8, v16, 5);
}


void PlaneWithCube::set_boundaries(Boundaries bounds) {
    m_bounds = bounds;
}


void PlaneWithCube::init_blocks() {
    check_params();

    double xl = 2 * m_xc - m_xmin;

    //снизу вверх
    v1 = BaseVertex::create(m_xmin, m_ymin, true);
    v2 = BaseVertex::create(m_xc,   m_ymin, false);
    v3 = BaseVertex::create(xl,     m_ymin, false);
    v4 = BaseVertex::create(m_xmax, m_ymin, true);

    v5 = BaseVertex::create(m_xc, m_yc - m_r, true);

    v6 = BaseVertex::create(m_xc - 0.5 * m_r, m_yc - 0.5 * m_r, false);
    v7 = BaseVertex::create(m_xc + 0.5 * m_r, m_yc - 0.5 * m_r, false);

    v8 = BaseVertex::create(m_xmin,      m_yc, false);
    v9 = BaseVertex::create(m_xc - m_r,  m_yc, true);
    v10 = BaseVertex::create(m_xc + m_r, m_yc, true);
    v11 = BaseVertex::create(xl,         m_yc, false);
    v12 = BaseVertex::create(m_xmax,     m_yc, false);

    v13 = BaseVertex::create(m_xc - 0.5 * m_r, m_yc + 0.5 * m_r, false);
    v14 = BaseVertex::create(m_xc + 0.5 * m_r, m_yc + 0.5 * m_r, false);

    v15 = BaseVertex::create(m_xc, m_yc + m_r, true);

    v16 = BaseVertex::create(m_xmin, m_ymax, true);
    v17 = BaseVertex::create(m_xc, m_ymax, false);
    v18 = BaseVertex::create(xl, m_ymax, false);
    v19 = BaseVertex::create(m_xmax, m_ymax, true);


    // Ограничивающие прямые области
    cube_side1 = Plane::create(v9, v5);
    cube_side2 = Plane::create(v5, v10);
    cube_side3 = Plane::create(v10, v15);
    cube_side4 = Plane::create(v15, v9);
    left   = Plane::create(v1, v16);
    right  = Plane::create(v4, v19);
    bottom = Plane::create(v1, v4);
    top    = Plane::create(v16, v19);

    cube_side1->set_boundary(m_bounds.hole);
    cube_side2->set_boundary(m_bounds.hole);
    cube_side3->set_boundary(m_bounds.hole);
    cube_side4->set_boundary(m_bounds.hole);
    left->set_boundary(m_bounds.left);
    right->set_boundary(m_bounds.right);
    bottom->set_boundary(m_bounds.bottom);
    top->set_boundary(m_bounds.top);

    // Генератор сетки
    m_blocks[0] = {v1, v6, v8, v9};
    m_blocks[0].set_boundary(v1, v8, left);
    m_blocks[0].set_boundary(v6, v9, cube_side1);

    m_blocks[1] = {v1, v2, v5, v6};
    m_blocks[1].set_boundary(v1, v2, bottom);
    m_blocks[1].set_boundary(v5, v6, cube_side1);

    m_blocks[2] = {v2, v3, v5, v7};
    m_blocks[2].set_boundary(v2, v3, bottom);
    m_blocks[2].set_boundary(v5, v7, cube_side2);

    m_blocks[3] = {v3, v7, v10, v11};
    m_blocks[3].set_boundary(v7, v10, cube_side2);

    m_blocks[4] = {v3, v4, v11, v12};
    m_blocks[4].set_boundary(v3, v4, bottom);
    m_blocks[4].set_boundary(v4, v12, right);

    m_blocks[5] = {v11, v12, v18, v19};
    m_blocks[5].set_boundary(v12, v19, right);
    m_blocks[5].set_boundary(v18, v19, top);

    m_blocks[6] = {v10, v11, v14, v18};
    m_blocks[6].set_boundary(v10, v14, cube_side3);

    m_blocks[7] = {v14, v15, v17, v18};
    m_blocks[7].set_boundary(v14, v15, cube_side3);
    m_blocks[7].set_boundary(v17, v18, top);

    m_blocks[8] = {v13, v15, v16, v17};
    m_blocks[8].set_boundary(v13, v15, cube_side4);
    m_blocks[8].set_boundary(v16, v17, top);

    m_blocks[9] = {v8, v9, v13, v16};
    m_blocks[9].set_boundary(v9, v13, cube_side4);
    m_blocks[9].set_boundary(v8, v16, left);

    // Необходимо связать блоки
    link();

    // Точность сглаживания (необязательно)
    set_accuracy(1.0e-5);
}

Box PlaneWithCube::bbox() const {
    Vector3d vmin(m_xmin, m_ymin, 0.0);
    Vector3d vmax(m_xmax, m_ymax, 0.0);

    return Box(vmin, vmax);
}

void PlaneWithCube::check_params() const {
    if (m_xmin >= m_xmax) {
        std::string message = "PlaneWithCube Error: x_min >= x_max";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (m_ymin >= m_ymax) {
        std::string message = "PlaneWithCube Error: y_min >= y_max";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }

    double R = m_r * m_xi;
    if (m_xmin >= m_xc - R || m_xmax <= m_xc + R
        || m_ymin >= m_yc - R || m_ymax <= m_yc + R) {
        std::string message = "PlaneWithCube Error: big hole";
        std::cerr << message << "\n";
        throw std::runtime_error(message);

    }
}

} // namespace zephyr::geom::generator::collection