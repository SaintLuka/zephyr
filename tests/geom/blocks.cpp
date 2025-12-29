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
    auto v1 = BaseVertex::create(-1.0, 0.0, true);
    auto v2 = BaseVertex::create(0.0, -0.5, true);
    auto v3 = BaseVertex::create(+1.0, 0.0, true);
    auto v4 = BaseVertex::create(0.0, +1.5, true);

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

// сектор с конформным модулем M
BlockStructured test2(double M) {
    double R = std::exp(M * M_PI / 2.0);

    auto v1 = BaseVertex::create(0.0, 1.0, true);
    auto v2 = BaseVertex::create(0.0, R, true);
    auto v3 = BaseVertex::create(1.0, 0.0, true);
    auto v4 = BaseVertex::create(R, 0.0, true);

    // Ограничивающие прямые области
    auto left  = Plane::create(v1, v2);
    auto right = Plane::create(v3, v4);
    //auto inner = Plane::create(v1, v3);
    //auto outer = Plane::create(v2, v4);
    auto inner = Circle::create(1.0, Vector3d::Zero());
    auto outer = Circle::create(R, Vector3d::Zero());

    left->set_boundary(Boundary::WALL);
    right->set_boundary(Boundary::ZOE);
    inner->set_boundary(Boundary::WALL);
    outer->set_boundary(Boundary::ZOE);

    BlockStructured blocks(4);

    // Генератор сетки
    blocks[0] = {v1, v2, v3, v4};
    blocks[0].set_boundary(v1, v2, left);
    blocks[0].set_boundary(v3, v4, right);
    blocks[0].set_boundary(v1, v3, inner);
    blocks[0].set_boundary(v2, v4, outer);

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

    auto v1 = BaseVertex::create(xmin, ymin, true);
    auto v3 = BaseVertex::create(xmax, ymin, true);
    auto v5 = BaseVertex::create(xmax, ymax, true);
    auto v7 = BaseVertex::create(xmin, ymax, true);

    auto v2 = BaseVertex::create(xc - r0 * std::sin(M_PI / 4.0),
                            yc - r0 * std::cos(M_PI / 4.0),
                            false);
    auto v4 = BaseVertex::create(xc + r0 * std::sin(M_PI / 4.0),
                            yc - r0 * std::cos(M_PI / 4.0),
                            false);
    auto v6 = BaseVertex::create(xc + r0 * std::sin(M_PI / 4.0),
                            yc + r0 * std::cos(M_PI / 4.0),
                            false);
    auto v8 = BaseVertex::create(xc - r0 * std::sin(M_PI / 4.0),
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
    auto v1 = BaseVertex::create(xmin, ymin, true);
    auto v2 = BaseVertex::create(x_l,    ymin, false);
    auto v3 = BaseVertex::create(x_r,    ymin, false);
    auto v4 = BaseVertex::create(xmax, ymin, true);

    auto v5 = BaseVertex::create(xmin, y_b, false);
    auto v6 = BaseVertex::create(x_l,    y_b, false);
    auto v7 = BaseVertex::create(x_r,    y_b, false);
    auto v8 = BaseVertex::create(xmax, y_b, false);

    auto v9  = BaseVertex::create(xc - a, yc - a, false);
    auto v10 = BaseVertex::create(xc + a, yc - a, false);
    auto v11 = BaseVertex::create(xc - a, yc + a, false);
    auto v12 = BaseVertex::create(xc + a, yc + a, false);

    auto v13 = BaseVertex::create(xmin, y_t, false);
    auto v14 = BaseVertex::create(x_l,    y_t, false);
    auto v15 = BaseVertex::create(x_r,    y_t, false);
    auto v16 = BaseVertex::create(xmax, y_t, false);

    auto v17 = BaseVertex::create(xmin, ymax, true);
    auto v18 = BaseVertex::create(x_l,    ymax, false);
    auto v19 = BaseVertex::create(x_r,    ymax, false);
    auto v20 = BaseVertex::create(xmax, ymax, true);

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
    BlockStructured blocks = test9();

    blocks.link();
    blocks.optimize();

    return 0;
}