#include "fast.h"

#include <zephyr/geom/generator/collection/wedge.h>
#include <zephyr/geom/generator/collection/semicircle_cutout.h>

#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/tests/shock-in-a-box.h>
#include <zephyr/phys/tests/mach.h>
#include <zephyr/phys/tests/sod.h>
#include <zephyr/phys/tests/blast_wave.h>
#include <zephyr/phys/tests/RiemannTest2D.h>
#include <zephyr/phys/tests/toro.h>

#include <zephyr/math/solver/riemann.h>
#include <zephyr/phys/eos/stiffened_gas.h>
#include <zephyr/math/solver/sm_fluid.h>

using zephyr::geom::generator::collection::Wedge;
using zephyr::geom::generator::collection::SemicircleCutout;
using zephyr::geom::generator::Rectangle;

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

using zephyr::math::RiemannSolver;
using zephyr::math::SmFluid;

// Для быстрого доступа по типу
SmFluid::State U;

/// Переменные для сохранения
double get_rho(AmrStorage::Item& cell) { return cell(U).rho; }
double get_u(AmrStorage::Item& cell) { return cell(U).v.x(); }
double get_v(AmrStorage::Item& cell) { return cell(U).v.y(); }
double get_w(AmrStorage::Item& cell) { return cell(U).v.z(); }
double get_p(AmrStorage::Item& cell) { return cell(U).p; }
double get_e(AmrStorage::Item& cell) { return cell(U).e; }


int main() {
    // Тестовая задача
    // Mach test(2.81);
    ToroTest test(100);
    // SodTest test;
    // BlastWave test(4.5);
    // RiemannTest2D test(6);

    // Уравнение состояния
    Eos& eos = test.eos;

    // Файл для записи
    PvdFile pvd("mesh", "/mnt/d/Quirk`s HLLC-LM"); //blastfail
    pvd.unique_nodes = true;

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"v", get_v};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"c",
                      [&eos](AmrStorage::Item& cell) -> double {
                          return eos.sound_speed_rp(cell(U).rho, cell(U).p);
                      }};

    

    Rectangle gen(0.0, 400.0, 0.0, 20.0);
    gen.set_nx(400);
    gen.set_ny(20);
    gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
                        .bottom=Boundary::WALL, .top=Boundary::WALL});

    // // Часть области с регулярной сеткой
    // auto fix_condition = [&test](const Vector3d& v) {
    //     // return v.x() <= test.x_jump + 0.01 * (test.xmax() - test.xmin());
    //     return v.x() <= 0.09;
    // };

    // Wedge gen(0.0, 3.0, 0.0, 1.0, 0.1666, 0.0, 
    //             {.left=Boundary::ZOE, .right=Boundary::ZOE,
    //             .bottom=Boundary::ZOE, .top=Boundary::ZOE});
    // gen.set_nx(40);
    // gen.set_fixed(fix_condition);


    //SemicircleCutout gen(0.49, 0.7, 0.0, 0.07, 0.6, 0.02);
    // SemicircleCutout gen(0.0, 0.21, 0.0, 0.07, 0.11, 0.02, 
    //                      {.left=Boundary::ZOE, .right=Boundary::ZOE, 
    //                       .bottom=Boundary::WALL, .top=Boundary::ZOE});
    // gen.set_ny(10);
    // gen.set_fixed(fix_condition);
    // gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
    //                     .bottom=Boundary::WALL, .top=Boundary::ZOE});

    // Создать сетку
    EuMesh mesh(U, &gen);
    int n_cells = mesh.n_cells();

    // Создать решатель
    auto solver = zephyr::math::SmFluid(eos, Fluxes::HLLC_LM);
    solver.set_accuracy(2);
    solver.set_CFL(0.4);

    // mesh.set_max_level(3);
    // mesh.set_distributor(solver.distributor());
    
    double time = 0.0;
    double next_write = 0.0;
    size_t n_step = 0;

    solver.init_cells(mesh, test);

    // for (int k = 0; k < mesh.max_level() + 3; ++k) {
    //     solver.init_cells(mesh, test);
    //     solver.set_flags(mesh);
    //     mesh.refine();
    // }

    while (time <= 1.01 * test.max_time()) {

        std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                  << "\tTime: " << std::setw(6) << std::setprecision(3) << time << "\n";
        if (time >= next_write) {
            pvd.save(mesh, time);
            next_write += test.max_time() / 100;
        };

        // Обновляем слои
        solver.update(mesh);

        // solver.set_flags(mesh);
        
        // mesh.refine();
        
        n_step += 1;
        time = solver.get_time();
    }

    auto fprint = [](const std::string &name, double value) {
        std::cout << name << ": " << value << '\n';
    };

    return 0;
}
