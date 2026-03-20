// Протестирую решатель на разных сетках

#include <iostream>
#include <iomanip>

#include <zephyr/phys/tests/test_2D.h>

#include <zephyr/math/solver/sm_fluid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/csv_file.h>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/threads.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/collection/wedge.h>
#include <zephyr/geom/generator/collection/plane_with_hole.h>

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;
using namespace zephyr::geom::generator;

using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::math::SmFluid;
using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;

auto WALL = Boundary::WALL;
auto ZOE = Boundary::ZOE;

// Простой квадрат с декартовой сеткой
Grid test1() {
    Rectangle gen(-1.0, 1.0, -1.0, 1.0);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_nx(200);
    return gen.make();
}

// Ячейки Вороного в прямоугольнике
Grid test2() {
    Rectangle gen(-1.0, 1.0, -1.0, 1.0, true);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_nx(200);
    return gen.make();
}

// Декартова сетка разбитая на треугольники (1)
Grid test3() {
    Rectangle gen(-1.0, 1.0, -1.0, 1.0);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_nx(100);
    Grid grid = gen.make();
    grid.triangulation(1);
    return grid;
}

// Декартова сетка разбитая на треугольники (2)
Grid test4() {
    Rectangle gen(-1.0, 1.0, -1.0, 1.0);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_nx(100);
    Grid grid = gen.make();
    grid.triangulation(2);
    return grid;
}

// Декартова сетка в кубе
Grid test5() {
    Cuboid gen(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL, .back=WALL, .front=WALL});
    gen.set_nx(40);
    return gen.make();
}

// Декартова сетка в кубе, побитая на пирамидки
Grid test6() {
    Cuboid gen(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL, .back=WALL, .front=WALL});
    gen.set_nx(20);
    Grid grid = gen.make();
    grid.pyramidize();
    return grid;
}

// Декартова сетка в кубе, побитая на тетраэдры
Grid test7() {
    Cuboid gen(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL, .back=WALL, .front=WALL});
    gen.set_nx(20);
    Grid grid = gen.make();
    grid.triangulation();
    return grid;
}

// Extrude декартовой сетки
Grid test8() {
    Rectangle gen(-1.0, 1.0, -1.0, 1.0);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_nx(40);
    Grid grid = gen.make();
    grid.extrude(Vector3d{0.1, 0.1, 1.0}, 20, WALL, WALL);
    return grid;
}

// Extrude сетки из треугольников
Grid test9() {
    Rectangle gen(-1.0, 1.0, -1.0, 1.0);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_nx(30);
    Grid grid = gen.make();
    grid.triangulation(2);
    grid.extrude(Vector3d{0.0, 0.0, 1.0}, 15, WALL, WALL);
    return grid;
}

// Extrude сетки из многоугольников
Grid test10() {
    Rectangle gen(-1.0, 1.0, -1.0, 1.0, true);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_nx(30);
    Grid grid = gen.make();
    grid.extrude(Vector3d{0.0, 0.0, 1.0}, 15, WALL, WALL);
    return grid;
}

// Простая BlockStructured сетка
Grid test11() {
    collection::PlaneWithHole gen(0.0, 2.0, 0.0, 1.0, 0.6, 0.2, 0.1);
    gen.set_boundaries({.left=ZOE, .right=ZOE, .bottom=WALL, .top=WALL});
    gen.set_ny(80);
    return gen.make();
}

// BlockStructured сетка с отражением и триангуляцией
Grid test12() {
    collection::Wedge gen(0.0, 1.0, 0.0, 0.5, 0.5, 0.2);
    gen.set_boundaries({.left=ZOE, .right=ZOE, .bottom=WALL, .top=WALL});
    gen.set_ny(20);

    Grid grid = gen.make();
    grid.triangulation(1);
    grid.mirror_x();
    grid.move({0.0, -0.5, 0.0});
    grid.mirror_y();
    return grid;
}

// BlockStructured с extrude
Grid test13() {
    collection::PlaneWithHole gen = collection::PlaneWithHole(0.0, 2.0, 0.0, 2.0, 1.0, 1.0, 0.3);
    gen.set_boundaries({.left=ZOE, .right=ZOE, .bottom=WALL, .top=WALL});
    gen.set_ny(80);
    Grid grid = gen.make();
    grid.extrude(Vector3d::UnitZ(), 40, ZOE, ZOE);
    return grid;
}

// BlockStructured, затем make_amr
Grid test14() {
    collection::PlaneWithHole gen = collection::PlaneWithHole(0.0, 2.0, 0.0, 2.0, 1.0, 1.0, 0.3);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_ny(80);
    Grid grid = gen.make();
    grid.make_amr();
    return grid;
}

// BlockStructured, затем extrude и make_amr
Grid test15() {
    collection::PlaneWithHole gen = collection::PlaneWithHole(0.0, 2.0, 0.0, 2.0, 1.0, 1.0, 0.3);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_ny(80);
    Grid grid = gen.make();
    grid.extrude(Vector3d::UnitZ()/20, 2, ZOE, ZOE);
    grid.make_amr();
    return grid;
}

int main(int argc, char** argv) {
    mpi::handler handler(argc, argv);
    threads::init(argc, argv);
    threads::info();
    threads::on();

    // Создать сетку
    Grid grid = test15();
    bool polyhedral = grid.polyhedral();
    EuMesh mesh(std::move(grid));

    // Создать и настроить решатель
    auto eos = IdealGas::create("Air");
    SmFluid solver(eos);
    solver.set_accuracy(2);
    solver.set_CFL(0.3);
    solver.set_limiter("MC");
    solver.set_method(Fluxes::HLLC);

    // Добавляем типы на сетку, выбираем основной слой
    auto data = solver.add_types(mesh);
    auto z = data.init;

    // Настройка сетки
    mesh.set_decomposition("XY");
    mesh.set_max_level(2);
    mesh.set_distributor(solver.distributor());

    // Задание начальных данных
    double R = 0.3;
    auto init_cells = [R, eos, z](EuMesh& mesh) {
        mesh.for_each([R, eos, z](EuCell& cell) {
            if (cell.center().norm() < R) {
                cell[z].density = 3.0;
                cell[z].pressure = 3.0;
            }
            else {
                cell[z].density = 1.0;
                cell[z].pressure = 1.0;
            }
            cell[z].velocity = Vector3d::Zero();
            cell[z].energy = eos->energy_rP(cell[z].density, cell[z].pressure);
        });
    };

    // Файл для записи
    PvdFile pvd("mesh", "output");
    //pvd.unique_nodes = true;
    pvd.polyhedral = polyhedral;

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables += {"density",  [z](EuCell& cell) -> double { return cell[z].density; }};
    pvd.variables += {"vel.x",    [z](EuCell& cell) -> double { return cell[z].velocity.x(); }};
    pvd.variables += {"vel.y",    [z](EuCell& cell) -> double { return cell[z].velocity.y(); }};
    pvd.variables += {"vel.z",    [z](EuCell& cell) -> double { return cell[z].velocity.z(); }};
    pvd.variables += {"pressure", [z](EuCell& cell) -> double { return cell[z].pressure; }};
    pvd.variables += {"energy",   [z](EuCell& cell) -> double { return cell[z].energy; }};

    double curr_time = 0.0;

    // Инициализация начальными данными
    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        init_cells(mesh);
        solver.set_flags(mesh);
        mesh.refine();
    }
    init_cells(mesh);

    size_t n_step = 0;
    double next_write = 0.0;

    Stopwatch elapsed(true);
    while (curr_time < 1.0) {
        if (curr_time >= next_write) {
            mpi::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(8) << std::setprecision(3) << curr_time << "\n";

            pvd.save(mesh, curr_time);
            next_write += 0.01;
        }

        // Обновляем слои
        solver.update(mesh);
        solver.set_flags(mesh);
        mesh.refine();

        curr_time += solver.dt();
        n_step += 1;
    }
    pvd.save(mesh, curr_time);
    elapsed.stop();

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    return 0;
}
