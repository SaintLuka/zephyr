#include <zephyr/geom/generator/curve/plane.h>
#include <zephyr/geom/generator/curve/circle.h>
#include <zephyr/geom/generator/block_structured.h>
#include <zephyr/geom/grid.h>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/vtu_file.h>

#include <zephyr/geom/generator/sector.h>
#include <zephyr/geom/generator/collection/wedge.h>
#include <zephyr/geom/generator/collection/semicircle_cutout.h>
#include <zephyr/geom/generator/collection/plane_with_hole.h>
#include <zephyr/geom/generator/collection/plane_with_cube.h>

#include "zephyr/geom/generator/rectangle.h"

using namespace zephyr::geom;
using namespace zephyr::mesh;
using namespace zephyr::io;
using namespace generator;

// Дельтоид, конформный модуль M = 1
BlockStructured test1() {
    auto v1 = BaseNode::create(-1.0, 0.0);
    auto v2 = BaseNode::create(0.0, -0.5);
    auto v3 = BaseNode::create(+1.0, 0.0);
    auto v4 = BaseNode::create(0.0, +1.5);

    // Генератор сетки
    BlockStructured blocks;

    blocks += {v1, v2, v3, v4};
    blocks.set_boundary(v1, v2, Boundary::WALL);
    blocks.set_boundary(v2, v3, Boundary::WALL);
    blocks.set_boundary(v3, v4, Boundary::WALL);
    blocks.set_boundary(v1, v4, Boundary::WALL);

    blocks.set_size(v1, v2, 20);
    return blocks;
}

// Сектор с конформным модулем M
BlockStructured test2(double M, double alpha=0.5*M_PI) {
    double R0 = 1.0;
    double R1 = std::exp(M * alpha);

    double sin = std::sin(0.5*alpha);
    double cos = std::cos(0.5*alpha);
    auto v1 = BaseNode::create(R0*cos, -R0*sin);
    auto v2 = BaseNode::create(R1*cos, -R1*sin);
    auto v3 = BaseNode::create(R0*cos, +R0*sin);
    auto v4 = BaseNode::create(R1*cos, +R1*sin);

    // Ограничивающие прямые области
    auto inner = Circle::create(1.0, Vector3d::Zero());
    auto outer = Circle::create(R1, Vector3d::Zero());

    // Генератор сетки
    BlockStructured blocks;

    blocks += {v1, v2, v3, v4};
    blocks.set_boundary(v1, v2, Boundary::WALL);
    blocks.set_boundary(v3, v4, Boundary::WALL);
    blocks.set_boundary(v1, v3, inner);
    blocks.set_boundary(v2, v4, outer);

    blocks.set_size(v1, v2, 50);
    return blocks;
}

// Сектор, разделенный на две части, конформные модули M1 и M2 у левой и
// правой частей. Позволяет протестировать сопряжение двух блоков.
// При зафиксированных точках на границах должен достигаться баланс.
BlockStructured test3(double M1, double M2, double angle=0.5*M_PI) {
    double R0 = 1.0;
    double R1 = std::exp(M1 * angle);
    double R2 = std::exp((M1 + M2) * angle);

    double sin = std::sin(0.5*angle);
    double cos = std::cos(0.5*angle);
    auto v1 = BaseNode::create(R0*cos, -R0*sin);
    auto v2 = BaseNode::create(R1*cos, -R1*sin, true);
    auto v3 = BaseNode::create(R2*cos, -R2*sin);
    auto v4 = BaseNode::create(R0*cos, +R0*sin);
    auto v5 = BaseNode::create(R1*cos, +R1*sin, true);
    auto v6 = BaseNode::create(R2*cos, +R2*sin);

    // Ограничивающие прямые области
    auto inner = Circle::create(R0, Vector3d::Zero());
    auto outer = Circle::create(R2, Vector3d::Zero());

    // Генератор сетки
    BlockStructured blocks;
    blocks += {v1, v2, v4, v5};
    blocks += {v2, v3, v5, v6};
    blocks.set_boundary({v4, v5, v6}, Boundary::WALL);
    blocks.set_boundary({v1, v2, v3}, Boundary::WALL);
    blocks.set_boundary(v1, v4, inner);
    blocks.set_boundary(v3, v6, outer);

    blocks.set_size(v1, v2, 30);
    return blocks;
}

// Сектор, разделенный на четыре части, у трёх частей заданы конформные
// модули M, Mx, My. Позволяет протестировать сопряжение четырёх блоков.
// При зафиксированных точках на границах должен достигаться баланс.
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

    auto v1 = BaseNode::create(R0*cos(phi1), R0*sin(phi1));
    auto v2 = BaseNode::create(R1*cos(phi1), R1*sin(phi1), true);
    auto v3 = BaseNode::create(R2*cos(phi1), R2*sin(phi1), true);
    auto v4 = BaseNode::create(R0*cos(phi2), R0*sin(phi2));
    auto v5 = BaseNode::create(R1*cos(phi2), R1*sin(phi2));
    auto v6 = BaseNode::create(R2*cos(phi2), R2*sin(phi2));
    auto v7 = BaseNode::create(R0*cos(phi3), R0*sin(phi3), true);
    auto v8 = BaseNode::create(R1*cos(phi3), R1*sin(phi3), true);
    auto v9 = BaseNode::create(R2*cos(phi3), R2*sin(phi3));

    // Ограничивающие прямые области
    auto inner = Circle::create(R0, Vector3d::Zero());
    auto outer = Circle::create(R2, Vector3d::Zero());

    // Генератор сетки
    BlockStructured blocks;
    blocks += {v1, v2, v4, v5};
    blocks += {v2, v3, v5, v6};
    blocks += {v4, v5, v7, v8};
    blocks += {v5, v6, v8, v9};
    blocks.set_boundary({v7, v8, v9}, Boundary::WALL);
    blocks.set_boundary({v1, v2, v3}, Boundary::WALL);
    blocks.set_boundary({v1, v4, v7}, inner);
    blocks.set_boundary({v3, v6, v9}, outer);

    blocks.set_size(v1, v2, 30);
    return blocks;
}

// Прямоугольник с круглым вырезом
BlockStructured test5() {
    double xmin =0.0, xmax = 2.5, ymin = -1.5, ymax = 1.5;
    double xc = 0.8, yc = 0.0, r0 = 0.2;

    double xi = 1.1;
    double x_l = xc - xi * r0;
    double x_r = xc + xi * r0;
    double y_b = yc - xi * r0;
    double y_t = yc + xi * r0;

    double a = r0 / std::sqrt(2.0);

    // Задаем базисные вершины для блоков
    auto v1 = BaseNode::create(xmin, ymin);
    auto v2 = BaseNode::create(x_l,  ymin);
    auto v3 = BaseNode::create(x_r,  ymin);
    auto v4 = BaseNode::create(xmax, ymin);

    auto v5 = BaseNode::create(xmin, y_b);
    auto v6 = BaseNode::create(x_l,  y_b);
    auto v7 = BaseNode::create(x_r,  y_b);
    auto v8 = BaseNode::create(xmax, y_b);

    auto v9  = BaseNode::create(xc - a, yc - a);
    auto v10 = BaseNode::create(xc + a, yc - a);
    auto v11 = BaseNode::create(xc - a, yc + a);
    auto v12 = BaseNode::create(xc + a, yc + a);

    auto v13 = BaseNode::create(xmin, y_t);
    auto v14 = BaseNode::create(x_l,  y_t);
    auto v15 = BaseNode::create(x_r,  y_t);
    auto v16 = BaseNode::create(xmax, y_t);

    auto v17 = BaseNode::create(xmin, ymax);
    auto v18 = BaseNode::create(x_l,  ymax);
    auto v19 = BaseNode::create(x_r,  ymax);
    auto v20 = BaseNode::create(xmax, ymax);

    // Генератор сетки
    BlockStructured blocks;
    blocks += {v1, v6, v2, v5};
    blocks += {v2, v3, v7, v6};
    blocks += {v3, v8, v4, v7};
    blocks += {v5, v13, v6, v14};
    blocks += {v9, v6, v11, v14};
    blocks += {v6, v7, v10, v9};
    blocks += {v7, v10, v15, v12};
    blocks += {v11, v15, v12, v14};
    blocks += {v7, v8, v16, v15};
    blocks += {v13, v18, v14, v17};
    blocks += {v14, v15, v19, v18};
    blocks += {v15, v16, v20, v19};

    blocks.set_boundary({v1, v5, v13, v17},   Boundary::ZOE);  // Левая граница
    blocks.set_boundary({v4, v8, v16, v20},   Boundary::ZOE);  // Правая граница
    blocks.set_boundary({v1, v2, v3, v4},     Boundary::WALL); // Нижняя граница
    blocks.set_boundary({v17, v18, v19, v20}, Boundary::WALL); // Верхняя граница

    // Окружность в центре области
    auto circle = Circle::create(r0, {xc, yc, 0.0});
    circle->set_boundary(Boundary::WALL);
    blocks.set_boundary({v9, v11, v12, v10, v9}, circle);

    // Установить по левой/правой границе
    blocks.set_size({v1, v5, v13, v17}, 160);
    blocks.set_size({v4, v8, v16, v20}, 160);

    // Установить по нижней/верхней границе
    //blocks.set_size({v1, v2, v3, v4}, 160);
    //blocks.set_size({v17, v18, v19, v20}, 160);

    // По радиусу кольца (допустим, более плоские ячейки)
    //blocks.set_size(v6, v9, 30);

    // Установить по внутреннему кольцу
    //blocks.set_size({v9, v10, v12, v11, v9}, 40);

    return blocks;
}

// Не связные блоки
BlockStructured test6() {
    double x1 = 0.0, x2 = 1.0, x3 = 2.0, x4 = 3.0;
    double y1 = 0.0, y2 = 1.0;

    // Задаем базисные вершины для структурированных блоков
    auto v1 = BaseNode::create(x1, y1);
    auto v2 = BaseNode::create(x2, y1);
    auto v3 = BaseNode::create(x2, y2);
    auto v4 = BaseNode::create(x1, y2);

    auto v5 = BaseNode::create(x3, y1);
    auto v6 = BaseNode::create(x4, y1);
    auto v7 = BaseNode::create(x4, y2);
    auto v8 = BaseNode::create(x3, y2);

    // Генератор сетки
    BlockStructured blocks;
    blocks += {v1, v2, v3, v4};
    blocks += {v5, v6, v7, v8};

    blocks.set_boundary({v1, v2, v3, v4, v1}, Boundary::WALL);
    blocks.set_boundary({v5, v6, v7, v8, v5}, Boundary::WALL);

    return blocks;
}

int main() {
    //BlockStructured blocks = test1();
    //BlockStructured blocks = test2(2.4, 0.3);
    //BlockStructured blocks = test3(0.5, 2.0, 0.3);
    //BlockStructured blocks = test4(1.0, 1.5, 0.25, 0.3);
    //BlockStructured blocks = test5();
    //BlockStructured blocks = test6();

    //blocks.set_verbosity(5);
    //blocks.set_iters_count(2000);
    //blocks.optimize();

    //Rectangle gen(0.0, 1.0, 0.0, 1.0);
    //gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE, .bottom=Boundary::WALL, .top=Boundary::WALL});
    //gen.set_ny(10);

    collection::Wedge gen(0.0, 2.0, 0.0, 1.0, 1.0, 0.3);
    gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE, .bottom=Boundary::WALL, .top=Boundary::WALL});
    gen.set_ny(20);

    //collection::SemicircleCutout gen(0.0, 2.5, 0.0, 1.2, 0.8, 0.2);
    //gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE, .bottom=Boundary::WALL, .top=Boundary::WALL});
    //gen.set_ny(80);

    //collection::PlaneWithHole gen(0.0, 2.0, 0.0, 1.0, 0.6, 0.2, 0.1);
    //gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE, .bottom=Boundary::WALL, .top=Boundary::WALL});
    //gen.set_ny(80);

    //collection::PlaneWithCube gen(0.0, 2.0, 0.0, 1.0, 1.0, 0.5, 0.1);
    //gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE, .bottom=Boundary::WALL, .top=Boundary::WALL});
    //gen.set_ny(80);

    Grid grid = gen.make();
    //grid.mirror_x();
    //grid.move(-Vector3d::UnitY());
    //grid.mirror_y();
    //grid.triangulation(2);
    //grid.make_amr();


    EuMesh mesh(std::move(grid));
    VtuFile::save("mesh.vtu", mesh, Variables{"faces2D"});

    return 0;
}