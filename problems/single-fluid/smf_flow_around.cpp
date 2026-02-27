/// @file smf_flow_around.cpp
/// @brief Двумерные газодинамические задачи с обтеканием сложных геометрий.
/// Или с взаимодействием со сложной геометрией.
#include <iostream>
#include <iomanip>

#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/collection/wedge.h>
#include <zephyr/geom/generator/collection/semicircle_cutout.h>
#include <zephyr/geom/generator/collection/plane_with_hole.h>
#include <zephyr/geom/generator/collection/plane_with_cube.h>
#include <zephyr/geom/generator/cuboid.h>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/math/solver/sm_fluid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/csv_file.h>

#include <zephyr/phys/tests/test_1D.h>

using zephyr::geom::generator::collection::Wedge;
using zephyr::geom::generator::collection::SemicircleCutout;
using zephyr::geom::generator::collection::PlaneWithHole;
using zephyr::geom::generator::collection::PlaneWithCube;
using zephyr::geom::generator::Rectangle;
using zephyr::geom::generator::Cuboid;

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

using zephyr::mesh::EuMesh;
using zephyr::math::SmFluid;
using zephyr::utils::threads;

int main() {
    threads::on();

    // Тестовая задача
    ShockWave test(3.0, 0.1, 1.0);

    // Уравнение состояния
    auto eos = test.get_eos();

    /*
    // Сеточные генераторы на выбор
    Rectangle gen(0.0, 1.0, 0.0, 0.2);
    gen.set_nx(200);
    gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
                        .bottom=Boundary::WALL, .top=Boundary::WALL});
    */

    /*
    Wedge gen(0.0, 3.0, 0.0, 1.8, 1.5, M_PI / 6.0);
    gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
                        .bottom=Boundary::WALL, .top=Boundary::WALL});
    gen.set_nx(100);
    */

    /*
    gen.set_fixed(fix_condition);

    // Часть области с регулярной сеткой
    auto fix_condition = [&test](const Vector3d& v) {
        return v.x() <= test.x_jump + 0.1 * (test.xmax() - test.xmin());
    };
     */


    // PlaneWithCube gen(0.0, 2.0, 0.0, 1.0, 0.6, 0.5, 0.1);
    // gen.set_boundaries({.left   = Boundary::ZOE, .right  = Boundary::ZOE,
    //                     .bottom = Boundary::ZOE, .top    = Boundary::ZOE,
    //                     .hole   = Boundary::WALL});
    // gen.set_nx(500);

    /*
    PlaneWithHole gen(0.0, 2.0, 0.0, 1.0, 0.6, 0.5, 0.1);
    gen.set_boundaries({.left   = Boundary::ZOE, .right  = Boundary::ZOE,
                        .bottom = Boundary::ZOE, .top    = Boundary::ZOE,
                        .hole   = Boundary::WALL});
    gen.set_nx(320);
    */

    Cuboid gen(0.0, 1.0,
               0.0, 1.0,
               0.0, 0.2);
    gen.set_boundaries({.left   = Boundary::ZOE, .right  = Boundary::ZOE,
                        .bottom = Boundary::ZOE, .top    = Boundary::ZOE,
                        .back   = Boundary::ZOE, .front  = Boundary::ZOE});
    gen.set_nx(20);
    gen.set_ny(20);
    gen.set_nz(20);


    //Wedge gen(0.0, 2.0, 0.0, 1.0, 1.0, 0.3);
    //gen.set_boundaries({.left=Boundary::WALL, .right=Boundary::WALL, .bottom=Boundary::WALL, .top=Boundary::WALL});
    //gen.set_ny(50);

    Grid grid = gen.make();
    //grid.mirror_x();
    //grid.move(-Vector3d::UnitY());
    //grid.mirror_y();
    //grid.triangulation(2);


    // Создать сетку
    EuMesh mesh(std::move(grid));

    // Создать решатель
    SmFluid solver(eos);
    solver.set_accuracy(2);
    solver.set_CFL(0.5);
    solver.set_limiter("MC");
    solver.set_method(Fluxes::HLLC_M);

    // Добавляем типы на сетку, выбираем основной слой
    auto data = solver.add_types(mesh);
    auto z = data.init;

    // Настройка сетки
    mesh.set_max_level(0);
    mesh.set_distributor(solver.distributor());

    // Начальные данные
    auto init_cell = [&test, z, eos](EuCell& cell) {
        auto cell_c = cell.center();
        if (cell_c.norm() < 0.2) {
            cell[z].density = 2.0;
            cell[z].pressure =3.0;
            cell[z].velocity = Vector3d::Zero();
        }
        else {
            cell[z].density = 1.0;
            cell[z].pressure =1.0;
            cell[z].velocity = Vector3d::Zero();
        }
        cell[z].energy = eos->energy_rP(cell[z].density, cell[z].pressure);
        //cell[z].density  = test.density(cell_c);
        //cell[z].velocity = test.velocity(cell_c);
        //cell[z].pressure = test.pressure(cell_c);
        //cell[z].energy   = test.energy(cell_c);
    };

    // Файл для записи
    PvdFile pvd("flow", "output");

    // Переменные для сохранения
    pvd.variables = {"level", "faces2D"};
    pvd.variables += {"density",  [z](EuCell& cell) -> double { return cell[z].density; }};
    pvd.variables += {"vel.x",    [z](EuCell& cell) -> double { return cell[z].velocity.x(); }};
    pvd.variables += {"vel.y",    [z](EuCell& cell) -> double { return cell[z].velocity.y(); }};
    pvd.variables += {"pressure", [z](EuCell& cell) -> double { return cell[z].pressure; }};
    pvd.variables += {"energy",   [z](EuCell& cell) -> double { return cell[z].energy; }};

    // Инициализация начальными данными
    if (mesh.adaptive()) {
        for (int k = 0; k < mesh.max_level() + 3; ++k) {
            mesh.for_each(init_cell);
            solver.set_flags(mesh);
            mesh.refine();
        }
    }
    mesh.for_each(init_cell);

    size_t n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    test.finish *= 15.0;

    while (curr_time < test.max_time()) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);
            next_write += test.max_time() / 100;
        }

        // Точное завершение в end_time
        solver.set_max_dt(test.max_time() - curr_time);

        // Обновляем слои
        solver.update(mesh);
        solver.set_flags(mesh);
        mesh.refine();

        curr_time += solver.dt();
        n_step += 1;
    }
    pvd.save(mesh, curr_time);

    return 0;
}