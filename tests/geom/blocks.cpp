#include <zephyr/geom/generator/curve/plane.h>
#include <zephyr/geom/generator/curve/circle.h>
#include <zephyr/geom/generator/block_structured.h>


using namespace zephyr::geom;
using namespace generator;

struct Boundaries {
    Boundary left, right, bottom, top, hole;
};

// Ромб, конформный модуль M = 1
BlockStructured test1() {
    auto v1 = BaseNode::create(-1.0, 0.0, true);
    auto v2 = BaseNode::create(0.0, -0.5, true);
    auto v3 = BaseNode::create(+1.0, 0.0, true);
    auto v4 = BaseNode::create(0.0, +1.5, true);

    // Ограничивающие прямые области
    auto left   = Plane::create(v1, v2);
    auto right  = Plane::create(v2, v3);
    auto bottom = Plane::create(v3, v4);
    auto top    = Plane::create(v1, v4);

    left->set_boundary(Boundary::WALL);
    right->set_boundary(Boundary::ZOE);
    bottom->set_boundary(Boundary::WALL);
    top->set_boundary(Boundary::ZOE);

    BlockStructured blocks(4);

    // Генератор сетки
    blocks[0] = {v1, v2, v3, v4};
    blocks[0].set_boundary(v1, v2, left);
    blocks[0].set_boundary(v2, v3, right);
    blocks[0].set_boundary(v3, v4, bottom);
    blocks[0].set_boundary(v1, v4, top);

    return blocks;
}

// Сектор с конформным модулем M
BlockStructured test2(double M, double alpha=0.5*M_PI) {
    double R0 = 1.0;
    double R1 = std::exp(M * alpha);

    double sin = std::sin(0.5*alpha);
    double cos = std::cos(0.5*alpha);
    auto v1 = BaseNode::create(R0*cos, -R0*sin, true);
    auto v2 = BaseNode::create(R1*cos, -R1*sin, true);
    auto v3 = BaseNode::create(R0*cos, +R0*sin, true);
    auto v4 = BaseNode::create(R1*cos, +R1*sin, true);

    // Ограничивающие прямые области
    auto right = Plane::create(v1, v2);
    auto left  = Plane::create(v3, v4);
    auto inner = Circle::create(1.0, Vector3d::Zero());
    auto outer = Circle::create(R1, Vector3d::Zero());

    left->set_boundary(Boundary::WALL);
    right->set_boundary(Boundary::ZOE);
    inner->set_boundary(Boundary::WALL);
    outer->set_boundary(Boundary::ZOE);

    BlockStructured blocks(4);

    // Генератор сетки
    blocks[0] = {v1, v2, v3, v4};
    blocks[0].set_boundary(v1, v2, right);
    blocks[0].set_boundary(v3, v4, left);
    blocks[0].set_boundary(v1, v3, inner);
    blocks[0].set_boundary(v2, v4, outer);

    return blocks;
}

// Сектор, разделенный на две части, конформные модули M1 и M2 у левой и
// правой частей. Позволяет протестировать сопряжение двух блоков.
BlockStructured test3(double M1, double M2, double angle=0.5*M_PI) {
    double R0 = 1.0;
    double R1 = std::exp(M1 * angle);
    double R2 = std::exp((M1 + M2) * angle);

    double sin = std::sin(0.5*angle);
    double cos = std::cos(0.5*angle);
    auto v1 = BaseNode::create(R0*cos, -R0*sin, true);
    auto v2 = BaseNode::create(R1*cos, -R1*sin, true);
    auto v3 = BaseNode::create(R2*cos, -R2*sin, true);
    auto v4 = BaseNode::create(R0*cos, +R0*sin, true);
    auto v5 = BaseNode::create(R1*cos, +R1*sin, true);
    auto v6 = BaseNode::create(R2*cos, +R2*sin, true);

    // Ограничивающие прямые области
    auto left  = Plane::create(v4, v6);
    auto right = Plane::create(v1, v3);
    auto inner = Circle::create(R0, Vector3d::Zero());
    auto outer = Circle::create(R2, Vector3d::Zero());

    left->set_boundary(Boundary::WALL);
    right->set_boundary(Boundary::ZOE);
    inner->set_boundary(Boundary::WALL);
    outer->set_boundary(Boundary::ZOE);

    BlockStructured blocks(2);

    // Генератор сетки
    blocks[0] = {v1, v2, v4, v5};
    blocks[0].set_boundary(v1, v4, inner);
    blocks[0].set_boundary(v4, v5, left);
    blocks[0].set_boundary(v1, v2, right);

    blocks[1] = {v2, v3, v5, v6};
    blocks[1].set_boundary(v2, v3, right);
    blocks[1].set_boundary(v5, v6, left);
    blocks[1].set_boundary(v3, v6, outer);

    return blocks;
}

// Сектор, разделенный на четыре части, у трёх частей заданы конформные
// модули M, Mx, My. Позволяет протестировать сопряжение четырёх блоков.
BlockStructured test4(double M, double Mx, double My, double angle=0.5*M_PI) {
    double x1 = 0.0, x2 = M, x3 = M + Mx;
    double y1 = 0.0, y2 = 1.0, y3 = (1.0 + My / M);

    // w(z) = exp((angle / y3) * z)
    // r_i = exp(x_i), phi_i = angle / y3 * y_i
    double R0 = std::exp((angle / y3) * x1);
    double R1 = std::exp((angle / y3) * x2);
    double R2 = std::exp((angle / y3) * x3);

    double phi1 = (angle / y3) * y1;
    double phi2 = (angle / y3) * y2;
    double phi3 = (angle / y3) * y3;

    using std::sin; using std::cos;

    auto v1 = BaseNode::create(R0*cos(phi1), R0*sin(phi1), true);
    auto v2 = BaseNode::create(R1*cos(phi1), R1*sin(phi1), true);
    auto v3 = BaseNode::create(R2*cos(phi1), R2*sin(phi1), true);
    auto v4 = BaseNode::create(R0*cos(phi2), R0*sin(phi2), true);
    auto v5 = BaseNode::create(R1*cos(phi2), R1*sin(phi2), false);
    auto v6 = BaseNode::create(R2*cos(phi2), R2*sin(phi2), true);
    auto v7 = BaseNode::create(R0*cos(phi3), R0*sin(phi3), true);
    auto v8 = BaseNode::create(R1*cos(phi3), R1*sin(phi3), true);
    auto v9 = BaseNode::create(R2*cos(phi3), R2*sin(phi3), true);

    // Ограничивающие прямые области
    auto left  = Plane::create(v7, v9);
    auto right = Plane::create(v1, v3);
    auto inner = Circle::create(R0, Vector3d::Zero());
    auto outer = Circle::create(R2, Vector3d::Zero());

    left->set_boundary(Boundary::WALL);
    right->set_boundary(Boundary::ZOE);
    inner->set_boundary(Boundary::WALL);
    outer->set_boundary(Boundary::ZOE);

    BlockStructured blocks(4);

    // Генератор сетки
    blocks[0] = {v1, v2, v4, v5};
    blocks[0].set_boundary(v1, v2, right);
    blocks[0].set_boundary(v1, v4, inner);

    blocks[1] = {v2, v3, v5, v6};
    blocks[1].set_boundary(v2, v3, right);
    blocks[1].set_boundary(v3, v6, outer);

    blocks[2] = {v4, v5, v7, v8};
    blocks[2].set_boundary(v4, v7, inner);
    blocks[2].set_boundary(v7, v8, left);

    blocks[3] = {v5, v6, v8, v9};
    blocks[3].set_boundary(v8, v9, left);
    blocks[3].set_boundary(v6, v9, outer);

    return blocks;
}

BlockStructured test8() {
    // params
    double xmin = 0.0, xmax = 3.0, ymin = 0.0, ymax = 1.0;
    double xc = 1.5, yc = 0.5, r0 = 0.3;

    Boundaries bounds = {
        Boundary::WALL, Boundary::WALL,
        Boundary::WALL, Boundary::WALL,
        Boundary::WALL
    };

    auto v1 = BaseNode::create(xmin, ymin, true);
    auto v3 = BaseNode::create(xmax, ymin, true);
    auto v5 = BaseNode::create(xmax, ymax, true);
    auto v7 = BaseNode::create(xmin, ymax, true);

    auto v2 = BaseNode::create(xc - r0 * std::sin(M_PI / 4.0),
                            yc - r0 * std::cos(M_PI / 4.0),
                            false);
    auto v4 = BaseNode::create(xc + r0 * std::sin(M_PI / 4.0),
                            yc - r0 * std::cos(M_PI / 4.0),
                            false);
    auto v6 = BaseNode::create(xc + r0 * std::sin(M_PI / 4.0),
                            yc + r0 * std::cos(M_PI / 4.0),
                            false);
    auto v8 = BaseNode::create(xc - r0 * std::sin(M_PI / 4.0),
                            yc + r0 * std::cos(M_PI / 4.0),
                            false);

    // Ограничивающие прямые области
    auto circle = Circle::create(r0, {xc, yc, 0.0});
    auto left   = Plane::create(v1, v7);
    auto right  = Plane::create(v3, v5);
    auto bottom = Plane::create(v1, v3);
    auto top    = Plane::create(v5, v7);

    left->set_boundary(bounds.left);
    right->set_boundary(bounds.right);
    bottom->set_boundary(bounds.bottom);
    top->set_boundary(bounds.top);
    circle->set_boundary(bounds.hole);

    BlockStructured blocks(4);

    // Генератор сетки
    blocks[0] = {v1, v2, v7, v8};
    blocks[0].set_boundary(v1, v7, left);
    blocks[0].set_boundary(v2, v8, circle);

    blocks[1] = {v1, v3, v2, v4};
    blocks[1].set_boundary(v1, v3, bottom);
    blocks[1].set_boundary(v2, v4, circle);

    blocks[2] = {v3, v4, v6, v5};
    blocks[2].set_boundary(v3, v5, right);
    blocks[2].set_boundary(v4, v6, circle);

    blocks[3] = {v5, v6, v7, v8};
    blocks[3].set_boundary(v5, v7, top);
    blocks[3].set_boundary(v6, v8, circle);


    return blocks;
}

BlockStructured test9() {
    double xmin =-1.6, xmax = 1.6, ymin = -0.6, ymax = +0.6;
    double xc = 0.0, yc = 0.0, r0 = 0.3;

    Boundaries bounds = {
        Boundary::WALL, Boundary::WALL,
        Boundary::WALL, Boundary::WALL,
        Boundary::WALL
    };
    

    double xi = 1.5;
    double x_l = xc - xi * r0;
    double x_r = xc + xi * r0;
    double y_b = yc - xi * r0;
    double y_t = yc + xi * r0;
    
    double a = r0 / std::sqrt(2.0);

    // Задаем базисные вершины для струтурированных блоков
    auto v1 = BaseNode::create(xmin, ymin, true);
    auto v2 = BaseNode::create(x_l,    ymin, false);
    auto v3 = BaseNode::create(x_r,    ymin, false);
    auto v4 = BaseNode::create(xmax, ymin, true);

    auto v5 = BaseNode::create(xmin, y_b, false);
    auto v6 = BaseNode::create(x_l,    y_b, false);
    auto v7 = BaseNode::create(x_r,    y_b, false);
    auto v8 = BaseNode::create(xmax, y_b, false);

    auto v9  = BaseNode::create(xc - a, yc - a, false);
    auto v10 = BaseNode::create(xc + a, yc - a, false);
    auto v11 = BaseNode::create(xc - a, yc + a, false);
    auto v12 = BaseNode::create(xc + a, yc + a, false);

    auto v13 = BaseNode::create(xmin, y_t, false);
    auto v14 = BaseNode::create(x_l,  y_t, false);
    auto v15 = BaseNode::create(x_r,  y_t, false);
    auto v16 = BaseNode::create(xmax, y_t, false);

    auto v17 = BaseNode::create(xmin, ymax, true);
    auto v18 = BaseNode::create(x_l,   ymax, false);
    auto v19 = BaseNode::create(x_r,    ymax, false);
    auto v20 = BaseNode::create(xmax, ymax, true);

    // Ограничивающие прямые области
    auto circle = Circle::create(r0, {xc, yc, 0.0});
    auto left   = Plane::create(v1, v17);
    auto right  = Plane::create(v4, v20);
    auto bottom = Plane::create(v1, v4);
    auto top    = Plane::create(v17, v20);

    left->set_boundary(bounds.left);
    right->set_boundary(bounds.right);
    bottom->set_boundary(bounds.bottom);
    top->set_boundary(bounds.top);
    circle->set_boundary(bounds.hole);
    
    
    BlockStructured blocks(12);

    // Генератор сетки
    blocks[0] = {v1, v2, v5, v6};
    blocks[0].set_boundary(v1, v5, left);
    blocks[0].set_boundary(v1, v2, bottom);

    blocks[1] = {v2, v3, v7, v6};
    blocks[1].set_boundary(v2, v3, bottom);

    blocks[2] = {v3, v4, v8, v7};
    blocks[2].set_boundary(v3, v4, bottom);
    blocks[2].set_boundary(v4, v8, right);

    blocks[3] = {v5, v6, v14, v13};
    blocks[3].set_boundary(v5, v13, left);

    blocks[4] = {v9, v6, v11, v14};
    blocks[4].set_boundary(v9, v11, circle);

    blocks[5] = {v6, v7, v10, v9};
    blocks[5].set_boundary(v9, v10, circle);

    blocks[6] = {v7, v10, v15, v12};
    blocks[6].set_boundary(v10, v12, circle);

    blocks[7] = {v11, v12, v15, v14};
    blocks[7].set_boundary(v11, v12, circle);

    blocks[8] = {v7, v8, v16, v15};
    blocks[8].set_boundary(v8, v16, right);

    blocks[9] = {v13, v14, v18, v17};
    blocks[9].set_boundary(v13, v17, left);
    blocks[9].set_boundary(v17, v18, top);

    blocks[10] = {v14, v15, v19, v18};
    blocks[10].set_boundary(v18, v19, top);

    blocks[11] = {v15, v16, v20, v19};
    blocks[11].set_boundary(v16, v20, right);
    blocks[11].set_boundary(v19, v20, top);
    
    /*
    // params

    auto v1 = BaseVertex::create(xc - r0 * std::sin(M_PI / 4.0), yc - r0 * std::cos(M_PI / 4.0), false);
    auto v2 = BaseVertex::create(xc + r0 * std::sin(M_PI / 4.0), yc - r0 * std::cos(M_PI / 4.0), false);
    auto v3 = BaseVertex::create(xc + r0 * std::sin(M_PI / 4.0), yc + r0 * std::cos(M_PI / 4.0), false);
    auto v4 = BaseVertex::create(xc - r0 * std::sin(M_PI / 4.0), yc + r0 * std::cos(M_PI / 4.0), false);

    auto v5 = BaseVertex::create(xmin, ymin, true);
    auto v6 = BaseVertex::create(x1, ymin, true);
    auto v7 = BaseVertex::create(x2, ymin, true);
    auto v8 = BaseVertex::create(xmax, ymin, true);

    // Ограничивающие прямые области
    auto circle = Circle::create(r0, {xc, yc, 0.0});
    auto left   = Plane::create(v1, v7);
    auto right  = Plane::create(v3, v5);
    auto bottom = Plane::create(v1, v3);
    auto top    = Plane::create(v5, v7);

    left->set_boundary(bounds.left);
    right->set_boundary(bounds.right);
    bottom->set_boundary(bounds.bottom);
    top->set_boundary(bounds.top);
    circle->set_boundary(bounds.hole);

    BlockStructured blocks(4);

    // Генератор сетки
    blocks[0] = {v1, v2, v7, v8};
    blocks[0].set_boundary(v1, v7, left);
    blocks[0].set_boundary(v2, v8, circle);

    blocks[1] = {v1, v3, v2, v4};
    blocks[1].set_boundary(v1, v3, bottom);
    blocks[1].set_boundary(v2, v4, circle);

    blocks[2] = {v3, v4, v6, v5};
    blocks[2].set_boundary(v3, v5, right);
    blocks[2].set_boundary(v4, v6, circle);

    blocks[3] = {v5, v6, v7, v8};
    blocks[3].set_boundary(v5, v7, top);
    blocks[3].set_boundary(v6, v8, circle);

    */
    return blocks;
}

int main() {
    //BlockStructured blocks = test1();
    //BlockStructured blocks = test2(2.4, 0.3);
    //BlockStructured blocks = test3(0.5, 2.0, 0.3);
    BlockStructured blocks = test4(1.0, 1.5, 0.25, 0.3);
    //BlockStructured blocks = test9();

    blocks.link();
    blocks.optimize();

    return 0;
}