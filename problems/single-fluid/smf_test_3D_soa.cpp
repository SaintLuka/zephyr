/// @file smf_test_3D.cpp
/// @brief Сферический взрыв на трехмерной декартовой сетке.

#include <iostream>
#include <iomanip>

#include <zephyr/mesh/mesh.h>
#include <zephyr/geom/generator/cuboid.h>

#include <zephyr/phys/tests/sedov_blast.h>

#include <zephyr/math/solver/sm_fluid_soa.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/threads.h>
#include <zephyr/utils/stopwatch.h>

#include "zephyr/mesh/euler/soa_mesh.h"
#include "zephyr/mesh/primitives/decomposition.h"

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

using zephyr::mesh::generator::Cuboid;
using zephyr::mesh::EuMesh;
using zephyr::mesh::QCell;
using zephyr::math::SmFluidSoA;

using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;
using zephyr::mesh::BFaces;

// Критерий адаптации подобран под задачу.
// Адаптация ячеек с плотностью выше 1.5.
void set_flags(SmFluidSoA& solver, SoaMesh &mesh) {
    if (!mesh.is_adaptive()) return;

    mesh.for_each([&solver](QCell& cell) {
        const double threshold = 1.5;

        if (cell(solver.curr).density > threshold) {
            cell.set_flag(1);
            return;
        }
        for (auto face: cell.faces()) {
            auto neib = face.neib();
            if (neib(solver.curr).density > threshold) {
                cell.set_flag(1);
                return;
            }
        }
        cell.set_flag(-1);
    });
}

int main(int argc, char** argv) {
    mpi::handler init(argc, argv);
    threads::init(argc, argv);
    threads::info();

    // Тестовая задача
    SedovBlast3D test({.gamma=1.4, .rho0=1.0, .E=1.0});

    // Уравнение состояния
    IdealGas::Ptr eos = test.eos;

    // Задание начальных данных
    double t0 = 0.1;   // test.time_by_radius(0.4);
    test.finish = 0.7; // test.time_by_radius(0.9);

    // Файл для записи
    PvdFile pvd("Sedov", "output");
    pvd.unique_nodes = false;

    size_t n_step = 0;
    double curr_time = t0;
    double next_write = t0;


    // Генератор сетки
    Cuboid gen = Cuboid(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    gen.set_nx(20);
    gen.set_boundaries({.left=Boundary::WALL, .right=Boundary::ZOE,
                        .bottom=Boundary::WALL, .top=Boundary::ZOE,
                        .back=Boundary::WALL, .front=Boundary::ZOE});

    SoaMesh mesh(gen);

    using zephyr::mesh::QCell;

    // Создать решатель
    SmFluidSoA solver(eos);
    solver.set_accuracy(2);
    solver.set_CFL(0.5);
    solver.set_limiter("MC");
    solver.set_method(Fluxes::HLLC_M);

    // Сеточная адаптация
    mesh.set_max_level(1);
    mesh.set_distributor(solver.distributor());

    // Переменные для сохранения
    pvd.variables = {"rank", "index", "level", "flag"};
    pvd.variables += {"rho", [&solver](QCell& cell) -> double { return cell(solver.curr).density; }};
    pvd.variables += {"v", [&solver](QCell& cell) -> double { return cell(solver.curr).velocity.norm(); }};
    pvd.variables += {"P", [&solver](QCell& cell) -> double { return cell(solver.curr).pressure; }};
    pvd.variables += {"e", [&solver](QCell& cell) -> double { return cell(solver.curr).energy; }};
    pvd.variables += {"rho.exact",
                      [&test, &curr_time](QCell &cell) -> double {
                          return test.density(cell.center().norm(), curr_time);
                      }};
    pvd.variables += {"v.exact",
                      [&test, &curr_time](QCell &cell) -> double {
                          return test.velocity(cell.center().norm(), curr_time);
                      }};
    pvd.variables += {"P.exact",
                      [&test, &curr_time](QCell &cell) -> double {
                          return test.pressure(cell.center().norm(), curr_time);
                      }};
    pvd.variables += {"e.exact",
                      [&test, &curr_time](QCell &cell) -> double {
                          return test.energy(cell.center().norm(), curr_time);
                      }};

    solver.curr = mesh.add_data<PState>("curr");
    solver.half = mesh.add_data<PState>("half");
    solver.next = mesh.add_data<PState>("next");
    solver.d_dx = mesh.add_data<PState>("d_dx");
    solver.d_dy = mesh.add_data<PState>("d_dy");
    solver.d_dz = mesh.add_data<PState>("d_dz");

    auto init_cells = [&](SoaMesh& mesh) {
        mesh.for_each([&](QCell& cell) {
            double r = cell.center().norm();
            Vector3d n = cell.center() / r;
            PState z;
            z.density  = test.density (r, t0);
            z.velocity = test.velocity(r, t0) * n;
            z.pressure = std::max(test.pressure(r, t0), 1.0e-3);
            z.energy   = eos->energy_rP(z.density, z.pressure);

            cell(solver.curr) = z;
        });
    };

    for (int k = 0; k < mesh.max_level() + 1; ++k) {
        pvd.save(mesh, k);
        init_cells(mesh);
        pvd.save(mesh, k+0.1);
        set_flags(solver, mesh);
        pvd.save(mesh, k + 0.2);
        mesh.refine();
        pvd.save(mesh, k + 0.3);
    }
    init_cells(mesh);
    pvd.save(mesh, 1000);
    return 0;


    Stopwatch sw_update;
    Stopwatch sw_flags;
    Stopwatch sw_refine;

    Stopwatch sw_step;
    Stopwatch sw_grad;
    Stopwatch sw_flux1;
    Stopwatch sw_flux2;

    pvd.save(mesh, 0.0);

    Stopwatch elapsed(true);
    while (n_step < 10000 && curr_time < test.max_time()) {
        if (n_step % 10 == 0) {
            mpi::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(10) << std::setprecision(5) << curr_time << "\n";
        }
        // Точное завершение в end_time
        solver.set_max_dt(test.max_time() - curr_time);

        // Обновляем слои
        sw_update.resume();

        sw_step.resume();
        solver.compute_dt(mesh);
        sw_step.stop();

        sw_grad.resume();
        solver.compute_grad(mesh);
        sw_grad.stop();

        sw_flux1.resume();
        solver.fluxes_stage1(mesh);
        sw_flux1.stop();

        sw_flux2.resume();
        solver.fluxes_stage2(mesh);
        solver.swap(mesh);
        sw_flux2.stop();

        sw_update.stop();

        sw_flags.resume();
        set_flags(solver, mesh);
        sw_flags.stop();

        // Для переноса по градиентам
        //sw_grad.resume();
        //solver.compute_grad(mesh);
        //sw_grad.stop();

        //sw_refine.resume();
        //mesh.refine();
        //sw_refine.stop();

        curr_time += solver.dt();
        n_step += 1;
    }
    elapsed.stop();
    pvd.save(mesh, 1.0);

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    mpi::cout << "  Update time:  " << sw_update.extended_time()
              << " ( " << sw_update.milliseconds() << " ms)\n";
    mpi::cout << "  Compute dt :  " << sw_step.extended_time()
              << " ( " << sw_step.milliseconds() << " ms)\n";
    mpi::cout << "  Gradient   :  " << sw_grad.extended_time()
              << " ( " << sw_grad.milliseconds() << " ms)\n";
    mpi::cout << "  Fluxes 1   :  " << sw_flux1.extended_time()
              << " ( " << sw_flux1.milliseconds() << " ms)\n";
    mpi::cout << "  Fluxes 2   :  " << sw_flux2.extended_time()
              << " ( " << sw_flux2.milliseconds() << " ms)\n";
    //mpi::cout << "  Flags  time:  " << sw_flags.extended_time()
    //          << " ( " << sw_flags.milliseconds() << " ms)\n";
    //mpi::cout << "  Gradient   :  " << sw_grad.extended_time()
    //          << " ( " << sw_grad.milliseconds() << " ms)\n";
    //mpi::cout << "  Refine time:  " << sw_refine.extended_time()
    //          << " ( " << sw_refine.milliseconds() << " ms)\n";

    return 0;
}
