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
        double xw, double phi, 
        Boundaries bounds) :
        BlockStructured(2),
        m_xmin(xmin), m_xmax(xmax),
        m_ymin(ymin), m_ymax(ymax),
        m_xw(xw), m_phi(phi),
        m_bounds(bounds) {

    m_name = "collection.wedge";
    init_blocks();
}

inline double sqr(double x) {
    return x * x;
}

void Wedge::set_nx(int Nx) {
    // Характерный размер ячейки
    double DX = (m_xmax - m_xmin) / Nx;

    // Положение клина на правой границе
    double m_yw = m_ymin + (m_xmax - m_xw) * std::tan(m_phi);

    // Площади блоков
    double S1 = (m_xw - m_xmin) * (m_ymax - m_ymin);
    double S2 = 0.5 * (m_xmax - m_xw) * (m_ymax - m_ymin + m_ymax - m_yw);

    // Характерное число ячеек области
    int N = std::max(2, int(std::ceil((S1 + S2) / sqr(DX))));
    int Ny = std::max(1, N / Nx);

    int Nx1 = std::max(1, int(std::round(Ny * (m_xw - m_xmin) / (m_ymax - m_ymin))));
    int Nx2 = std::max(1, Nx - Nx1);

    // Нет необходимости устанавливать все размеры
    // у каждого блока, поскольку они связаны
    m_blocks[0].set_size(v1, v4, Ny);
    m_blocks[0].set_size(v1, v2, Nx1);
    m_blocks[1].set_size(v2, v3, Nx2);
}


void Wedge::set_boundaries(Boundaries bounds) {
    m_bounds = bounds;
}

void Wedge::init_blocks() {
    check_params();

    double m_yw = m_ymin + (m_xmax - m_xw) * std::tan(m_phi);

    // Задаем базисные вершины для струтурированных блоков
    v1 = BaseVertex::create(m_xmin, m_ymin, true);
    v2 = BaseVertex::create(m_xw,   m_ymin, true);
    v3 = BaseVertex::create(m_xmax, m_yw,   true);
    v4 = BaseVertex::create(m_xmin, m_ymax, true);
    v5 = BaseVertex::create(m_xw,   m_ymax, false);
    v6 = BaseVertex::create(m_xmax, m_ymax, true);

    // Ограничивающие прямые области
    left   = Plane::create(v1, v4);
    right  = Plane::create(v3, v6);
    bottom = Plane::create(v1, v2);
    top    = Plane::create(v4, v6);
    wedge  = Plane::create(v2, v3);

    left->set_boundary(m_bounds.left);
    right->set_boundary(m_bounds.right);
    bottom->set_boundary(m_bounds.bottom);
    top->set_boundary(m_bounds.top);
    // wedge->set_boundary(m_bounds.bottom);
    wedge->set_boundary(Boundary::WALL);

    // Генератор сетки
    m_blocks[0] = {v1, v2, v4, v5};
    m_blocks[0].set_boundary(v1, v4, left);
    m_blocks[0].set_boundary(v1, v2, bottom);
    m_blocks[0].set_boundary(v4, v5, top);

    m_blocks[1] = {v2, v3, v5, v6};
    m_blocks[1].set_boundary(v2, v3, wedge);
    m_blocks[1].set_boundary(v5, v6, top);
    m_blocks[1].set_boundary(v3, v6, right);

    // Необходимо связать блоки
    link();

    // Точность сглаживания (необязательно)
    set_accuracy(1.0e-5);
}


Box Wedge::bbox() const {   
    Vector3d vmin(m_xmin, m_ymin, 0.0);
    Vector3d vmax(m_xmax, m_ymax, 0.0);

    return Box(vmin, vmax);
}

void Wedge::check_params() const {
    if (m_xmin >= m_xmax) {
        std::string message = "Wedge Error: x_min >= x_max";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (m_ymin >= m_ymax) {
        std::string message = "Wedge Error: y_min >= y_max";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (m_xw >= m_xmax || m_xw <= m_xmin) {
        std::string message = "Wedge Error: x_w not in [x_min, x_max]";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (m_phi < 0.0 || m_phi > 0.75 * M_PI) {
        std::string message = "Wedge Error: phi not in [0.0, pi/2]";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }

    // высота клина
    double H = (m_xmax - m_xw) * std::tan(m_phi);
    if (H > 0.5 * (m_ymax - m_ymin)) {
        std::string message = "Wedge Error: big wedge, increase (y_max - y_min)";
        std::cerr << message << "\n";
        throw std::runtime_error(message);

    }
}

} // namespace zephyr::geom::generator::collection