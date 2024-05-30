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

/// @brief PlaneWithCube gen(0, 3.2, 0, 3.2, 1.6, 1.6, 0.3); должно быть симметричным. пока только так)
PlaneWithCube::PlaneWithCube(
        double xmin, double xmax,
        double ymin, double ymax,
        double xc, double yc,
        double r) :
        BlockStructured(8),
        m_xmin(xmin), m_xmax(xmax),
        m_ymin(ymin), m_ymax(ymax),
        m_xc(xc), m_yc(yc),
        m_r(r), m_xi(2.0) {

    m_name = "collection.plane-with-cube";

    init_blocks();
}

/// @brief 8 x 20 
void PlaneWithCube::set_nx(int Nx) {

    // TODO
    m_blocks[0].set_size(v1, v7, 4);
    m_blocks[0].set_size(v1, v5, 10);

    m_blocks[1].set_size(v1, v2, 4);

    m_blocks[2].set_size(v2, v3, 4);

    m_blocks[3].set_size(v3, v10, 4);

    m_blocks[4].set_size(v10, v16, 4);

    m_blocks[5].set_size(v15, v16, 4);

    m_blocks[6].set_size(v14, v15, 4);

    m_blocks[7].set_size(v7, v14, 4);
}


void PlaneWithCube::set_boundaries(Boundaries bounds) {
    m_bounds = bounds;
}


void PlaneWithCube::init_blocks() {
    check_params();

    //снизу вверх
    v1 = BaseVertex::create(m_xmin, m_ymin, true);
    v2 = BaseVertex::create(m_xc,   m_ymin, false);
    v3 = BaseVertex::create(m_xmax, m_ymin, true);

    v4 = BaseVertex::create(m_xc, m_yc - m_r, true);

    v5 = BaseVertex::create(m_xc - 0.5 * m_r, m_yc - 0.5 * m_r, false);
    v6 = BaseVertex::create(m_xc + 0.5 * m_r, m_yc - 0.5 * m_r, false);

    v7 = BaseVertex::create(m_xmin,      m_yc, false);
    v8 = BaseVertex::create(m_xc - m_r,  m_yc, true);
    v9 = BaseVertex::create(m_xc + m_r,  m_yc, true);
    v10 = BaseVertex::create(m_xmax,     m_yc, false);

    v11 = BaseVertex::create(m_xc - 0.5 * m_r, m_yc + 0.5 * m_r, false);
    v12 = BaseVertex::create(m_xc + 0.5 * m_r, m_yc + 0.5 * m_r, false);

    v13 = BaseVertex::create(m_xc, m_yc + m_r, true);

    v14 = BaseVertex::create(m_xmin, m_ymax, true);
    v15 = BaseVertex::create(m_xc,   m_ymax, false);
    v16 = BaseVertex::create(m_xmax, m_ymax, true);


    // Ограничивающие прямые области
    cube_side1 = Plane::create(v4, v8);
    cube_side2 = Plane::create(v4, v9);
    cube_side3 = Plane::create(v9, v13);
    cube_side4 = Plane::create(v13, v8);
    left   = Plane::create(v1, v14);
    right  = Plane::create(v3, v16);
    bottom = Plane::create(v1, v3);
    top    = Plane::create(v14, v16);

    cube_side1->set_boundary(m_bounds.hole);
    cube_side2->set_boundary(m_bounds.hole);
    cube_side3->set_boundary(m_bounds.hole);
    cube_side4->set_boundary(m_bounds.hole);
    left->set_boundary(m_bounds.left);
    right->set_boundary(m_bounds.right);
    bottom->set_boundary(m_bounds.bottom);
    top->set_boundary(m_bounds.top);

    // Генератор сетки
    m_blocks[0] = {v1, v5, v7, v8};
    m_blocks[0].set_boundary(v1, v7, left);
    m_blocks[0].set_boundary(v5, v8, cube_side1);

    m_blocks[1] = {v1, v2, v4, v5};
    m_blocks[1].set_boundary(v1, v2, bottom);
    m_blocks[1].set_boundary(v4, v5, cube_side1);

    m_blocks[2] = {v2, v3, v4, v6};
    m_blocks[2].set_boundary(v2, v3, bottom);
    m_blocks[2].set_boundary(v4, v6, cube_side2);

    m_blocks[3] = {v3, v6, v9, v10};
    m_blocks[3].set_boundary(v6, v9, cube_side2);
    m_blocks[3].set_boundary(v3, v10, right);

    m_blocks[4] = {v9, v10, v12, v16};
    m_blocks[4].set_boundary(v9, v12, cube_side3);
    m_blocks[4].set_boundary(v10, v16, right);

    m_blocks[5] = {v12, v13, v15, v16};
    m_blocks[5].set_boundary(v12, v13, cube_side3);
    m_blocks[5].set_boundary(v15, v16, top);

    m_blocks[6] = {v11, v13, v14, v15};
    m_blocks[6].set_boundary(v11, v13, cube_side4);
    m_blocks[6].set_boundary(v14, v15, top);

    m_blocks[7] = {v7, v8, v11, v14};
    m_blocks[7].set_boundary(v8, v11, cube_side4);
    m_blocks[7].set_boundary(v7, v14, left);

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