#include <iostream>
#include <algorithm>
#include <cassert>
#include <memory>

#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>
#include <zephyr/utils/json.h>

#include <zephyr/geom/generator/curve/plane.h>
#include <zephyr/geom/generator/curve/circle.h>
#include <zephyr/geom/generator/collection/plane_with_hole.h>

namespace zephyr::geom::generator::collection {

PlaneWithHole::PlaneWithHole(const Json& config) :
        Generator("collection.plane-with-hole"),
        m_xmin(0.0), m_xmax(0.0),
        m_ymin(0.0), m_ymax(0.0),
        m_xc(0.0), m_yc(0.0),
        m_r(0.0) {

    throw std::runtime_error("PlaneWithHole(const Json& config): not implemented");
/*
    if (!config["geometry"]) {
        throw std::runtime_error("EuMesh config doesn't contain 'geometry'");
    }

    m_xmin = config["geometry"]["x_min"].as<double>();
    m_xmax = config["geometry"]["x_max"].as<double>();
    m_ymin = config["geometry"]["y_min"].as<double>();
    m_ymax = config["geometry"]["y_max"].as<double>();
    m_xc   = config["geometry"]["x_c"].as<double>();
    m_yc   = config["geometry"]["y_c"].as<double>();
    m_r    = config["geometry"]["r"].as<double>();

    if (!config["boundary"]) {
        throw std::runtime_error("EuMesh config doesn't contain 'boundary'");
    }
    m_left_flag   = boundary_from_string(config["boundary"]["left"].as<std::string>());
    m_right_flag  = boundary_from_string(config["boundary"]["right"].as<std::string>());
    m_bottom_flag = boundary_from_string(config["boundary"]["bottom"].as<std::string>());
    m_top_flag    = boundary_from_string(config["boundary"]["top"].as<std::string>());
    m_hole_flag   = boundary_from_string(config["boundary"]["hole"].as<std::string>());

    init_blocks();

    int N = 0;
    if (config["cells_per_x"]) {
        N = config["cells_per_x"].as<int>();
    }
    else if (config["nx"]) {
        N = config["nx"].as<int>();
    }
    else {
        throw std::runtime_error("Config.mesh doesn't contain size 'nx'");
    }
    set_nx(N);

    if (config["verbose"]) {
        m_blocks.set_verbosity(config["verbosity"].as<bool>());
    }
    if (config["iters_count"]) {
        m_blocks.set_iters_count(config["iters_count"].as<int>());
    }
    */
}

PlaneWithHole::PlaneWithHole(
        double xmin, double xmax,
        double ymin, double ymax,
        double xc, double yc,
        double r) :
        Generator("collection.plane-with-hole"),
        m_xmin(xmin), m_xmax(xmax),
        m_ymin(ymin), m_ymax(ymax),
        m_xc(xc), m_yc(yc),
        m_r(r) {

    double xi = 1.1;
    double x_l = m_xc - xi * m_r;
    double x_r = m_xc + xi * m_r;
    double y_b = m_yc - xi * m_r;
    double y_t = m_yc + xi * m_r;

    double a = m_r / std::sqrt(2.0);

    // Задаем базисные вершины для блоков
    v1 = BaseNode::create(m_xmin, m_ymin);
    v2 = BaseNode::create(x_l,  m_ymin);
    v3 = BaseNode::create(x_r,  m_ymin);
    v4 = BaseNode::create(m_xmax, m_ymin);

    v5 = BaseNode::create(m_xmin, y_b);
    v6 = BaseNode::create(x_l,  y_b);
    v7 = BaseNode::create(x_r,  y_b);
    v8 = BaseNode::create(m_xmax, y_b);

    v9  = BaseNode::create(m_xc - a, m_yc - a);
    v10 = BaseNode::create(m_xc + a, m_yc - a);
    v11 = BaseNode::create(m_xc - a, m_yc + a);
    v12 = BaseNode::create(m_xc + a, m_yc + a);

    v13 = BaseNode::create(m_xmin, y_t);
    v14 = BaseNode::create(x_l,  y_t);
    v15 = BaseNode::create(x_r,  y_t);
    v16 = BaseNode::create(m_xmax, y_t);

    v17 = BaseNode::create(m_xmin, m_ymax);
    v18 = BaseNode::create(x_l,  m_ymax);
    v19 = BaseNode::create(x_r,  m_ymax);
    v20 = BaseNode::create(m_xmax, m_ymax);

    // Генератор сетки
    m_blocks += {v1, v6, v2, v5};
    m_blocks += {v2, v3, v7, v6};
    m_blocks += {v3, v8, v4, v7};
    m_blocks += {v5, v13, v6, v14};
    m_blocks += {v9, v6, v11, v14};
    m_blocks += {v6, v7, v10, v9};
    m_blocks += {v7, v10, v15, v12};
    m_blocks += {v11, v15, v12, v14};
    m_blocks += {v7, v8, v16, v15};
    m_blocks += {v13, v18, v14, v17};
    m_blocks += {v14, v15, v19, v18};
    m_blocks += {v15, v16, v20, v19};

    // Окружность в центре области
    circle = Circle::create(m_r, {m_xc, m_yc, 0.0});
    left   = Plane::create(v1, v17);
    right  = Plane::create(v4, v20);
    bottom = Plane::create(v1, v4);
    top    = Plane::create(v17, v20);

    m_blocks.set_boundary({v9, v11, v12, v10, v9}, circle);
    m_blocks.set_boundary({v1, v5, v13, v17}, left);
    m_blocks.set_boundary({v4, v8, v16, v20}, right);
    m_blocks.set_boundary({v1, v2, v3, v4}, bottom);
    m_blocks.set_boundary({v17, v18, v19, v20}, top);

    //m_blocks.plot_layout();
    //m_blocks.set_iters_count(1);
    m_blocks.optimize();
    //m_blocks.plot();
}

void PlaneWithHole::set_nx(int Nx) {
    // Установить по нижней/верхней границе
    m_blocks.set_size({v1, v2, v3, v4}, Nx);
    m_blocks.set_size({v17, v18, v19, v20}, Nx);
}

void PlaneWithHole::set_ny(int Ny) {
    // Установить по левой/правой границе
    m_blocks.set_size({v1, v5, v13, v17}, Ny);
    m_blocks.set_size({v4, v8, v16, v20}, Ny);
}

void PlaneWithHole::set_boundaries(Boundaries bounds) const {
    circle->set_boundary(bounds.hole);
    left->set_boundary(bounds.left);
    right->set_boundary(bounds.right);
    top->set_boundary(bounds.top);
    bottom->set_boundary(bounds.bottom);
}

Box PlaneWithHole::bbox() const {
    Vector3d vmin(m_xmin, m_ymin, 0.0);
    Vector3d vmax(m_xmax, m_ymax, 0.0);

    return {vmin, vmax};
}

void PlaneWithHole::check_params() const {
    if (m_xmin >= m_xmax) {
        throw std::runtime_error("PlaneWithHole Error: x_min >= x_max");
    }
    if (m_ymin >= m_ymax) {
        throw std::runtime_error("PlaneWithHole Error: y_min >= y_max");
    }

    double R = 1.2 * m_r;
    if (m_xmin >= m_xc - R || m_xmax <= m_xc + R || m_ymin >= m_yc - R || m_ymax <= m_yc + R) {
        throw std::runtime_error("PlaneWithHole Error: big hole");
    }
}

} // namespace zephyr::geom::generator::collection