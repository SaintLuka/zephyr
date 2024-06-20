#include "fast.h"

#include <zephyr/geom/generator/collection/wedge.h>
#include <zephyr/geom/generator/collection/semicircle_cutout.h>
#include <zephyr/geom/generator/collection/plane_with_hole.h>
#include <zephyr/geom/generator/collection/plane_with_cube.h>
#include <zephyr/geom/generator/cuboid.h>

#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/tests/mach.h>
#include <zephyr/phys/tests/sod.h>
#include <zephyr/phys/tests/sedov.h>
#include <zephyr/phys/tests/RiemannTest2D.h>
#include <zephyr/phys/tests/toro.h>
#include <zephyr/phys/tests/supersonic_flow_around_cylinder.h>


#include <zephyr/math/solver/riemann.h>
#include <zephyr/phys/eos/stiffened_gas.h>
#include <zephyr/math/solver/sm_fluid.h>

using zephyr::geom::generator::collection::Wedge;
using zephyr::geom::generator::collection::SemicircleCutout;
using zephyr::geom::generator::collection::PlaneWithHole;
using zephyr::geom::generator::collection::PlaneWithCube;
using zephyr::geom::generator::Rectangle;
using zephyr::geom::generator::Cuboid;

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

    threads::on(16);

    // Тестовая задача
    Mach test(2.85, 1.275, 0.5);
    // ToroTest test(1);
    // SodTest test;
    // RiemannTest2D test(6);
    // SuperSonicFlowAroundCylinder test(3);
    // SedovBlast test(4.5, 0.5);

    // Уравнение состояния
    Eos& eos = test.eos;

    // Файл для записи
    PvdFile pvd("mesh", "/mnt/c/cube_plane2"); //blastfail
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
    
    // Rectangle gen(test.xmin(), test.xmax(), test.ymin(), test.ymax());
    // gen.set_nx(4);
    // gen.set_ny(4);
    // gen.set_boundaries({.left=Boundary::WALL, .right=Boundary::WALL,
    //                     .bottom=Boundary::WALL, .top=Boundary::WALL});

    // Часть области с регулярной сеткой
    auto fix_condition = [&test](const Vector3d& v) {
        return v.x() <= test.x_jump + 0.01 * (test.xmax() - test.xmin());
    };

    //PlaneWithCube gen(0, 2.2, 0, 3.2, 0.7, 1.6, 0.3); // в статье // 1.375
    PlaneWithCube gen(0, 3.2, 0, 3.2, 1.6, 1.6, 0.3);
    gen.set_nx(-1);

    // Wedge gen(0.0, 3.0, 0.0, 1.0, 0.1666, 0.0, 
    //             {.left=Boundary::ZOE, .right=Boundary::ZOE,
    //             .bottom=Boundary::ZOE, .top=Boundary::ZOE});
    // gen.set_nx(40);
    // gen.set_fixed(fix_condition);

    // PlaneWithHole gen(test.xmin(), test.xmax(), test.ymin(), test.ymax(), 
    //                   0.3, 0.5 * (test.ymin() + test.ymax()), 0.1);
    // gen.set_boundaries({.left   = Boundary::ZOE, .right  = Boundary::ZOE,    
    //                     .bottom = Boundary::ZOE, .top    = Boundary::ZOE,
    //                     .hole   = Boundary::WALL});
    // gen.set_nx(20);

    // Cuboid gen(test.xmin(), test.xmax(), 
    //            test.ymin(), test.ymax(), 
    //            test.zmin(), test.zmax());
    // gen.set_boundaries({.left   = Boundary::ZOE, .right  = Boundary::ZOE,
    //                     .bottom = Boundary::ZOE, .top    = Boundary::ZOE,
    //                     .back   = Boundary::ZOE, .front  = Boundary::ZOE});
    // gen.set_nx(30);
    // gen.set_ny(30);
    // gen.set_nz(30);

    // Создать сетку
    EuMesh mesh(U, &gen);
    int n_cells = mesh.n_cells();

    // Создать решатель
    auto solver = zephyr::math::SmFluid(eos, Fluxes::HLLC);
    solver.set_accuracy(2);
    solver.set_CFL(0.4);

    mesh.set_max_level(5);
    mesh.set_distributor(solver.distributor());
    
    double time = 0.0;
    double next_write = 0.0;
    size_t n_step = 0;

    // solver.init_cells(mesh, test);

    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        solver.init_cells(mesh, test);
        solver.set_flags(mesh);
        mesh.refine();
    }

    while (time <= 1.01 * test.max_time()) {

        std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                  << "\tTime: " << std::setw(6) << std::setprecision(3) << time << "\n";
        if (time >= next_write) {
            pvd.save(mesh, time);
            next_write += test.max_time() / 100;
        };

        // Обновляем слои
        solver.update(mesh);
        solver.set_flags(mesh);
        mesh.refine();
        
        n_step += 1;
        time = solver.get_time();
    }

    auto fprint = [](const std::string &name, double value) {
        std::cout << name << ": " << value << '\n';
    };

    threads::off();

    return 0;
}
