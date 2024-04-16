// #include <iostream>
// #include <algorithm>
// #include <cassert>
// #include <memory>

// #include <zephyr/geom/box.h>
// #include <zephyr/geom/grid.h>

// #include <zephyr/geom/generator/collection/wedge2.h>
// #include <zephyr/geom/generator/curve/plane.h>
// #include <zephyr/geom/generator/curve/curve.h>

// namespace zephyr::geom::generator::collection {

// Wedge2::Wedge2(
//         double xmin, double xmax,
//         double ymin, double ymax,
//         double xc, double phi) :
//         BlockStructured(5),
//         m_xmin(xmin), m_xmax(xmax),
//         m_ymin(ymin), m_ymax(ymax),
//         m_xc(xc), m_phi(phi),
//         m_bounds() {

//     m_name = "collection.wedge2";

//     init_blocks();
// }

// void Wedge2::set_ny(int Ny) {
//     // Характерный размер ячейки
//     double H = (m_ymax - m_ymin);
//     double h = H / Ny;

//     // Характерный размер ячейки
//     double dx = (m_xmax - m_xmin) / Nx;

//     // Положение клина на правой границе
//     double m_yw = m_ymin + (m_xmax - m_xw) * std::tan(m_phi);

//     // // Нет необходимости устанавливать все размеры
//     // // у каждого блока, поскольку они связаны
//     // m_blocks[0].set_size(v1, v9, Ny);
//     // m_blocks[0].set_size(v1, v2, Nxl);

//     // m_blocks[2].set_size(v4, v10, Nr);
//     // m_blocks[2].set_size(v4, v5, Nxc);

//     // m_blocks[4].set_size(v7, v8, Nxr);
//     // m_blocks[4].set_size(v8, v12, Ny);

//     // Нет необходимости устанавливать все размеры
//     // у каждого блока, поскольку они связаны
//     m_blocks[0].set_size(v1, v4, Nxc);
//     m_blocks[0].set_size(v1, v2, Nxl);
//     m_blocks[1].set_size(v2, v3, Nxr);
// }


// void Wedge2::set_boundaries(Boundaries bounds) {
//     m_bounds = bounds;
// }

// void Wedge2::init_blocks() {
//     check_params();

//     double H = m_ymax - m_ymin;
//     double L =  m_xmax - m_xmin;

//     double x1 = xmax - x_mc;
//     double y1 = x1 * std::tan(phi);
//     double x2 = ymax - x_mc; 

//     // Ограничивающие прямые области
//     left   = Plane::create(v1, v4);
//     right  = Plane::create(v3, v6);
//     bottom = Plane::create(v1, v2);
//     top    = Plane::create(v4, v6);
//     wedge  = Plane::create(v2, v3);

//     left->set_boundary(m_bounds.left);
//     right->set_boundary(m_bounds.right);
//     bottom->set_boundary(m_bounds.bottom);
//     top->set_boundary(m_bounds.top);
//     wedge->set_boundary(m_bounds.bottom);

//     // Задаем базисные вершины для струтурированных блоков
//     v1 = BaseVertex::create(m_xmin, m_ymin, true);
//     v2 = BaseVertex::create(m_xc,   m_ymin, true);
//     v3 = BaseVertex::create(m_xmax, m_w,   true);
//     v4 = BaseVertex::create(m_xmin, m_ymax, true);
//     v5 = BaseVertex::create(m_xw,   m_ymax, false);
//     v6 = BaseVertex::create(m_xmax, m_ymax, true);

//     // Ограничивающие прямые области
//     left   = Plane::create(v1, v4);
//     right  = Plane::create(v3, v6);
//     bottom = Plane::create(v1, v2);
//     top    = Plane::create(v4, v6);
//     wedge  = Plane::create(v2, v3);

//     left->set_boundary(m_bounds.left);
//     right->set_boundary(m_bounds.right);
//     bottom->set_boundary(m_bounds.bottom);
//     top->set_boundary(m_bounds.top);
//     wedge->set_boundary(m_bounds.bottom);

//     // Генератор сетки
//     m_blocks[0] = {v1, v2, v4, v5};
//     m_blocks[0].set_boundary(v1, v4, left);
//     m_blocks[0].set_boundary(v1, v2, bottom);
//     m_blocks[0].set_boundary(v4, v5, top);

//     m_blocks[1] = {v2, v3, v5, v6};
//     m_blocks[1].set_boundary(v2, v3, wedge);
//     m_blocks[1].set_boundary(v5, v6, top);
//     m_blocks[1].set_boundary(v3, v6, right);

//     // Необходимо связать блоки
//     link();

//     // Точность сглаживания (необязательно)
//     set_accuracy(1.0e-5);
// }


// Box Wedge2::bbox() const {
//     Vector3d vmin(m_xmin, m_ymin, 0.0);
//     Vector3d vmax(m_xmax, m_ymax, 0.0);

//     return Box(vmin, vmax);
// }

// void Wedge2::check_params() const {
//     if (m_xmin >= m_xmax) {
//         std::string message = "Wedge2 Error: x_min >= x_max";
//         std::cerr << message << "\n";
//         throw std::runtime_error(message);
//     }
//     if (m_ymin >= m_ymax) {
//         std::string message = "Wedge2 Error: y_min >= y_max";
//         std::cerr << message << "\n";
//         throw std::runtime_error(message);
//     }
//     if (m_xc >= m_xmax || m_xc <= m_xmin) {
//         std::string message = "Wedge2 Error: x_w not in [x_min, x_max]";
//         std::cerr << message << "\n";
//         throw std::runtime_error(message);
//     }
//     if (2.99 * m_r > (m_ymax - m_ymin)) {
//         // По оси y умещается 3 радиуса
//         std::string message = "Wedge2 Error: big radius, increase (y_max - y_min)";
//         std::cerr << message << "\n";
//         throw std::runtime_error(message);
//     }
//     if (2.99 * m_r > (m_xmax - m_xc)) {
//         // По оси x умещается 3 радиуса вправо
//         std::string message = "Wedge2 Error: big radius, increase (x_max - x_c)";
//         std::cerr << message << "\n";
//         throw std::runtime_error(message);
//     }
//     if (2.99 * m_r > (m_xc - m_xmin)) {
//         // По оси x умещается 3 радиуса влево
//         std::string message = "Wedge2 Error: big radius, increase (x_c - x_min)";
//         std::cerr << message << "\n";
//         throw std::runtime_error(message);
//     }

//     if (m_ymax - m_ymin > (m_xmax - m_xc)) {
//         // По оси x умещается толщина пластины вправо
//         std::string message = "Wedge2 Error: small rectangle, increase (x_max - x_c)";
//         std::cerr << message << "\n";
//         throw std::runtime_error(message);
//     }
//     if (m_ymax - m_ymin > (m_xc - m_xmin)) {
//         // По оси x умещается толщина пластины влево
//         std::string message = "Wedge2 Error: small rectangle, increase (x_c - x_min)";
//         std::cerr << message << "\n";
//         throw std::runtime_error(message);
//     }
// }

// } // namespace zephyr::geom::generator::collection