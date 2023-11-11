#include <iostream>
#include <algorithm>
#include <cassert>
#include <memory>

#include <zephyr/geom/geom.h>
#include <zephyr/geom/box.h>
#include <zephyr/geom/grid.h>

#include <zephyr/geom/generator/curve/plane.h>
#include <zephyr/geom/generator/curve/circle.h>
#include <zephyr/geom/generator/collection/plane_with_hole.h>

namespace zephyr::geom::generator::collection {

#ifdef ZEPHYR_ENABLE_YAML
PlaneWithHole::PlaneWithHole(YAML::Node config) :
        Generator("collection.plane-with-hole"),
        blocks(12),
        m_xmin(0.0), m_xmax(0.0),
        m_ymin(0.0), m_ymax(0.0),
        m_xc(0.0), m_yc(0.0),
        m_r(0.0), m_xi(2.0),
        m_left_flag(Boundary::UNDEFINED), m_right_flag(Boundary::UNDEFINED),
        m_bottom_flag(Boundary::UNDEFINED), m_top_flag(Boundary::UNDEFINED),
        m_hole_flag(Boundary::UNDEFINED) {

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

    if (config["geometry"]["xi"]) {
        m_xi = config["geometry"]["xi"].as<double>();
    }

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
        std::string message = "Config.mesh doesn't contain size 'nx'";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    set_nx(N);

    if (config["verbose"]) {
        blocks.set_verbose(config["verbose"].as<bool>());
    }
    if (config["epsilon"]) {
        blocks.set_accuracy(config["epsilon"].as<double>());
    }
}
#endif

PlaneWithHole::PlaneWithHole(
        double xmin, double xmax,
        double ymin, double ymax,
        double xc, double yc,
        double r) :
        BlockStructured(12),
        m_xmin(xmin), m_xmax(xmax),
        m_ymin(ymin), m_ymax(ymax),
        m_xc(xc), m_yc(yc),
        m_r(r), m_xi(2.0),
        m_bounds() {

    m_name = "collection.plane-with-hole";

    init_blocks();
}

void PlaneWithHole::set_nx(int Nx) {
    // Характерный размер ячейки
    double DX = (m_xmax - m_xmin) / Nx;

    int Ny = static_cast<int>((m_ymax - m_ymin) / DX);
    int np = static_cast<int>(M_PI_2 * m_r * m_xi / DX);

    np = std::max(np, 1);
    Nx = std::max(Nx, np + 2);
    Ny = std::max(Ny, np + 2);

    int nr = static_cast<int>(std::log(m_xi) / std::log(1.0 + M_PI_2/np));

    int Nx1, Nx2;

    double dx1 = m_xc - m_xmin;
    double dx2 = m_xmax - m_xc;
    if (dx1 > dx2) {
        Nx1 = static_cast<int>(std::round(dx1/(dx1 + dx2) * Nx - np / 2.0));
        Nx1 = std::max(Nx1, 1);
        Nx2 = Nx - np - Nx1;
    }
    else {
        Nx2 = static_cast<int>(std::round(dx2/(dx1 + dx2) * Nx - np / 2.0));
        Nx2 = std::max(Nx2, 1);
        Nx1 = Nx - np - Nx2;
    }

    int Ny1, Ny2;

    double dy1 = m_yc - m_ymin;
    double dy2 = m_ymax - m_yc;
    if (dy1 > dy2) {
        Ny1 = static_cast<int>(std::round(dy1/(dy1 + dy2) * Ny - np / 2.0));
        Ny1 = std::max(Ny1, 1);
        Ny2 = Ny - np - Ny1;
    }
    else {
        Ny2 = static_cast<int>(std::round(dy2/(dy1 + dy2) * Ny - np / 2.0));
        Ny2 = std::max(Ny2, 1);
        Ny1 = Ny - np - Ny2;
    }

    if (Nx1 + np + Nx2 != Nx) {
        throw std::runtime_error("PlaneWithHole error: Strange error #1");
    }
    if (Ny1 + np + Ny2 != Ny) {
        throw std::runtime_error("PlaneWithHole error: Strange error #1");
    }

    // Нет необходимости устанавливать все размеры
    // у каждого блока, поскольку они связаны
    m_blocks[0].set_size(v1, v2, Nx1);
    m_blocks[1].set_size(v2, v3, np);
    m_blocks[2].set_size(v3, v4, Nx2);
    m_blocks[10].set_size(v18, v19, np);

    m_blocks[0].set_size(v1, v5, Ny1);
    m_blocks[3].set_size(v5, v13, np);
    m_blocks[9].set_size(v13, v17, Ny2);
    m_blocks[8].set_size(v8, v16, np);

    m_blocks[4].set_size(v6, v9, nr);
}


void PlaneWithHole::set_boundaries(Boundaries bounds) {
    m_bounds = bounds;
}

void PlaneWithHole::init_blocks() {
    check_params();

    double xi = 1.5;
    double x_l = m_xc - xi * m_r;
    double x_r = m_xc + xi * m_r;
    double y_b = m_yc - xi * m_r;
    double y_t = m_yc + xi * m_r;

    double a = m_r / std::sqrt(2.0);

    // Задаем базисные вершины для струтурированных блоков
    v1 = BaseVertex::create(m_xmin, m_ymin, true);
    v2 = BaseVertex::create(x_l,    m_ymin, false);
    v3 = BaseVertex::create(x_r,    m_ymin, false);
    v4 = BaseVertex::create(m_xmax, m_ymin, true);

    v5 = BaseVertex::create(m_xmin, y_b, false);
    v6 = BaseVertex::create(x_l,    y_b, false);
    v7 = BaseVertex::create(x_r,    y_b, false);
    v8 = BaseVertex::create(m_xmax, y_b, false);

    v9  = BaseVertex::create(m_xc - a, m_yc - a, false);
    v10 = BaseVertex::create(m_xc + a, m_yc - a, false);
    v11 = BaseVertex::create(m_xc - a, m_yc + a, false);
    v12 = BaseVertex::create(m_xc + a, m_yc + a, false);

    v13 = BaseVertex::create(m_xmin, y_t, false);
    v14 = BaseVertex::create(x_l,    y_t, false);
    v15 = BaseVertex::create(x_r,    y_t, false);
    v16 = BaseVertex::create(m_xmax, y_t, false);

    v17 = BaseVertex::create(m_xmin, m_ymax, true);
    v18 = BaseVertex::create(x_l,    m_ymax, false);
    v19 = BaseVertex::create(x_r,    m_ymax, false);
    v20 = BaseVertex::create(m_xmax, m_ymax, true);

    // Ограничивающие прямые области
    circle = Circle::create(m_r, {m_xc, m_yc, 0.0});
    left   = Plane::create(v1, v17);
    right  = Plane::create(v4, v20);
    bottom = Plane::create(v1, v4);
    top    = Plane::create(v17, v20);

    left->set_boundary(m_bounds.left);
    right->set_boundary(m_bounds.right);
    bottom->set_boundary(m_bounds.bottom);
    top->set_boundary(m_bounds.top);
    circle->set_boundary(m_bounds.hole);

    // Генератор сетки
    m_blocks[0] = {v1, v2, v5, v6};
    m_blocks[0].set_boundary(v1, v5, left);
    m_blocks[0].set_boundary(v1, v2, bottom);

    m_blocks[1] = {v2, v3, v7, v6};
    m_blocks[1].set_boundary(v2, v3, bottom);

    m_blocks[2] = {v3, v4, v8, v7};
    m_blocks[2].set_boundary(v3, v4, bottom);
    m_blocks[2].set_boundary(v4, v8, right);

    m_blocks[3] = {v5, v6, v14, v13};
    m_blocks[3].set_boundary(v5, v13, left);

    m_blocks[4] = {v6, v9, v11, v14};
    m_blocks[4].set_boundary(v9, v11, circle);

    m_blocks[5] = {v6, v7, v10, v9};
    m_blocks[5].set_boundary(v9, v10, circle);

    m_blocks[6] = {v10, v7, v15, v12};
    m_blocks[6].set_boundary(v10, v12, circle);

    m_blocks[7] = {v11, v12, v15, v14};
    m_blocks[7].set_boundary(v11, v12, circle);

    m_blocks[8] = {v7, v8, v16, v15};
    m_blocks[8].set_boundary(v8, v16, right);

    m_blocks[9] = {v13, v14, v18, v17};
    m_blocks[9].set_boundary(v13, v17, left);
    m_blocks[9].set_boundary(v17, v18, top);

    m_blocks[10] = {v14, v15, v19, v18};
    m_blocks[10].set_boundary(v18, v19, top);

    m_blocks[11] = {v15, v16, v20, v19};
    m_blocks[11].set_boundary(v16, v20, right);
    m_blocks[11].set_boundary(v19, v20, top);

    // Необходимо связать блоки
    link();

    // Точность сглаживания (необязательно)
    set_accuracy(1.0e-5);
}


Box PlaneWithHole::bbox() const {
    Vector3d vmin(m_xmin, m_ymin, 0.0);
    Vector3d vmax(m_xmax, m_ymax, 0.0);

    return Box(vmin, vmax);
}

void PlaneWithHole::check_params() const {
    if (m_xmin >= m_xmax) {
        std::string message = "PlaneWithHole Error: x_min >= x_max";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (m_ymin >= m_ymax) {
        std::string message = "PlaneWithHole Error: y_min >= y_max";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }

    double R = m_r * m_xi;
    if (m_xmin >= m_xc - R || m_xmax <= m_xc + R
        || m_ymin >= m_yc - R || m_ymax <= m_yc + R) {
        std::string message = "PlaneWithHole Error: big hole";
        std::cerr << message << "\n";
        throw std::runtime_error(message);

    }
}

} // namespace zephyr::geom::generator::collection