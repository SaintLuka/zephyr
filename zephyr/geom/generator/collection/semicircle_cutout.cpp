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
        double xc, double r, Boundaries bounds) :
        BlockStructured(5),
        m_xmin(xmin), m_xmax(xmax),
        m_ymin(ymin), m_ymax(ymax),
        m_xc(xc), m_r(r),
        m_bounds(bounds) {

    m_name = "collection.semicircle-cutout";

    init_blocks();
}

inline double sqr(double x) {
    return x * x;
}

void SemicircleCutout::set_ny(int Ny) {
    // Характерный размер ячейки
    double H = (m_ymax - m_ymin);
    double h = H / Ny;

    double xi = H / m_r; // Соотношение "окружностей"
    int Nr = static_cast<int>(std::log(xi) / std::log(1.0 + M_PI_4 / Ny));

    int Nxc = 2 * Ny;
    int Nxl = std::max(1, int(((m_xc - H) - m_xmin) / h));
    int Nxr = std::max(1, int((m_xmax - (m_xc + H)) / h));

    // Нет необходимости устанавливать все размеры
    // у каждого блока, поскольку они связаны
    m_blocks[0].set_size(v1, v9, Ny);
    m_blocks[0].set_size(v1, v2, Nxl);

    m_blocks[2].set_size(v4, v10, Nr);
    m_blocks[2].set_size(v4, v5, Nxc);

    m_blocks[4].set_size(v7, v8, Nxr);
    m_blocks[4].set_size(v8, v12, Ny);
}


void SemicircleCutout::set_boundaries(Boundaries bounds) {
    m_bounds = bounds;
}

void SemicircleCutout::init_blocks() {
    check_params();

    double H = m_ymax - m_ymin;

    double x1 = m_xc - H;
    double x6 = m_xc + H;
    assert(m_xmin < x1 && x6 < m_xmax);

    double x2 = m_xc - m_r;
    double x5 = m_xc + m_r;

    double x3 = m_xc - m_r * std::cos(M_PI / 3.0);
    double x4 = m_xc + m_r * std::cos(M_PI / 3.0);
    double yc = m_ymin + m_r * std::sin(M_PI / 3.0);
    assert(yc < m_ymax);


    // Задаем базисные вершины для структурированных блоков
    v1 = BaseVertex::create(m_xmin, m_ymin, true);
    v2 = BaseVertex::create(x1,     m_ymin, false);
    v3 = BaseVertex::create(x2,     m_ymin, true);
    v4 = BaseVertex::create(x3,     yc,     false);
    v5 = BaseVertex::create(x4,     yc,     false);
    v6 = BaseVertex::create(x5,     m_ymin, true);
    v7 = BaseVertex::create(x6,     m_ymin, false);
    v8 = BaseVertex::create(m_xmax, m_ymin, true);
    v9 = BaseVertex::create(m_xmin,  m_ymax, true);
    v10 = BaseVertex::create(x1,     m_ymax, false);
    v11 = BaseVertex::create(x6,     m_ymax, false);
    v12 = BaseVertex::create(m_xmax, m_ymax, true);

    // Ограничивающие прямые области
    left   = Plane::create(v1, v9);
    right  = Plane::create(v8, v12);
    bottom = Plane::create(v1, v8);
    top    = Plane::create(v9, v12);
    circle = Circle::create(v3, v4, v6);

    left->set_boundary(m_bounds.left);
    right->set_boundary(m_bounds.right);
    bottom->set_boundary(m_bounds.bottom);
    top->set_boundary(m_bounds.top);
    circle->set_boundary(m_bounds.bottom);

    // Генератор сетки
    m_blocks[0] = {v1, v2, v9, v10};
    m_blocks[0].set_boundary(v1, v9, left);
    m_blocks[0].set_boundary(v1, v2, bottom);
    m_blocks[0].set_boundary(v9, v10, top);

    m_blocks[1] = {v2, v3, v4, v10};
    m_blocks[1].set_boundary(v2, v3, bottom);
    m_blocks[1].set_boundary(v3, v4, circle);

    m_blocks[2] = {v4, v5, v10, v11};
    m_blocks[2].set_boundary(v4, v5, circle);
    m_blocks[2].set_boundary(v10, v11, top);

    m_blocks[3] = {v6, v7, v5, v11};
    m_blocks[3].set_boundary(v6, v7, bottom);
    m_blocks[3].set_boundary(v5, v6, circle);

    m_blocks[4] = {v7, v8, v11, v12};
    m_blocks[4].set_boundary(v7, v8, bottom);
    m_blocks[4].set_boundary(v11, v12, top);
    m_blocks[4].set_boundary(v8, v12, right);

    // Необходимо связать блоки
    link();

    // Точность сглаживания (необязательно)
    set_accuracy(1.0e-5);
}


Box SemicircleCutout::bbox() const {
    Vector3d vmin(m_xmin, m_ymin, 0.0);
    Vector3d vmax(m_xmax, m_ymax, 0.0);

    return Box(vmin, vmax);
}

void SemicircleCutout::check_params() const {
    if (m_xmin >= m_xmax) {
        std::string message = "SemicircleCutout Error: x_min >= x_max";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (m_ymin >= m_ymax) {
        std::string message = "SemicircleCutout Error: y_min >= y_max";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (m_xc >= m_xmax || m_xc <= m_xmin) {
        std::string message = "SemicircleCutout Error: x_w not in [x_min, x_max]";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (2.99 * m_r > (m_ymax - m_ymin)) {
        // По оси y умещается 3 радиуса
        std::string message = "SemicircleCutout Error: big radius, increase (y_max - y_min)";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (2.99 * m_r > (m_xmax - m_xc)) {
        // По оси x умещается 3 радиуса вправо
        std::string message = "SemicircleCutout Error: big radius, increase (x_max - x_c)";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (2.99 * m_r > (m_xc - m_xmin)) {
        // По оси x умещается 3 радиуса влево
        std::string message = "SemicircleCutout Error: big radius, increase (x_c - x_min)";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }

    if (m_ymax - m_ymin > (m_xmax - m_xc)) {
        // По оси x умещается толщина пластины вправо
        std::string message = "SemicircleCutout Error: small rectangle, increase (x_max - x_c)";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
    if (m_ymax - m_ymin > (m_xc - m_xmin)) {
        // По оси x умещается толщина пластины влево
        std::string message = "SemicircleCutout Error: small rectangle, increase (x_c - x_min)";
        std::cerr << message << "\n";
        throw std::runtime_error(message);
    }
}

} // namespace zephyr::geom::generator::collection