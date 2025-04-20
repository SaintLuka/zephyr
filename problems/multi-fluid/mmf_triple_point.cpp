/// @file Двумерная задача переноса

#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/strip.h>
#include <zephyr/mesh/mesh.h>

#include <zephyr/phys/tests/test_1D.h>

#include <zephyr/math/solver/riemann.h>
#include <zephyr/math/solver/mm_fluid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/csv_file.h>

#include <zephyr/utils/stopwatch.h>
#include <zephyr/utils/mpi.h>

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::mmf;

using zephyr::mesh::EuMesh;

using zephyr::utils::Stopwatch;
using zephyr::utils::threads;
using zephyr::utils::mpi;


// Для быстрого доступа по типу
MmFluid::State U;

/// Переменные для сохранения
double get_rho(AmrStorage::Item& cell) { return cell(U).density; }
double get_vx(AmrStorage::Item& cell) { return cell(U).velocity.x(); }
double get_vy(AmrStorage::Item& cell) { return cell(U).velocity.y(); }
double get_p(AmrStorage::Item& cell) { return cell(U).pressure; }
double get_e(AmrStorage::Item& cell) { return cell(U).energy; }
double get_T(AmrStorage::Item &cell) { return cell(U).temperature; }
double get_cln(AmrStorage::Item &cell) { return cell(U).mass_frac.index(); }
double get_mfrac1(AmrStorage::Item &cell) { return cell(U).mass_frac[0]; }
double get_mfrac2(AmrStorage::Item &cell) { return cell(U).mass_frac[1]; }
double get_vfrac1(AmrStorage::Item &cell) { return cell(U).vol_frac(0); }
double get_vfrac2(AmrStorage::Item &cell) { return cell(U).vol_frac(1); }
double get_rho1(AmrStorage::Item &cell) { return cell(U).densities[0]; }
double get_rho2(AmrStorage::Item &cell) { return cell(U).densities[1]; }
double get_normal_x(AmrStorage::Item &cell) { return cell(U).n[0].x(); }
double get_normal_y(AmrStorage::Item &cell) { return cell(U).n[0].y(); }


int main(int argc, char** argv) {
    mpi::handler init(argc, argv);
    threads::off();

    // Материал
    Eos::Ptr sg1 = IdealGas::create(1.4, 1.0);
    Eos::Ptr sg2 = IdealGas::create(1.5, 1.0);
    Eos::Ptr sg3 = IdealGas::create(1.5, 1.0);

    // Формальная смесь
    MixturePT mixture = {sg1, sg2, sg3};

    // Файл для записи
    PvdFile pvd("TP", "output");
    PvdFile pvd_body0("body0", "output");
    PvdFile pvd_body1("body1", "output");
    PvdFile pvd_body2("body2", "output");

    // Переменные для сохранения
    pvd.variables += {"cln", get_cln};
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"vx", get_vx};
    pvd.variables += {"vy", get_vy};
    pvd.variables += {"e", get_e};
    pvd.variables += {"P", get_p};
    pvd.variables += {"T", get_T};
    pvd.variables += {"b1", get_mfrac1};
    pvd.variables += {"b2", get_mfrac2};
    pvd.variables += {"a1", get_vfrac1};
    pvd.variables += {"a2", get_vfrac2};
    pvd.variables += {"rho1", get_rho1};
    pvd.variables += {"rho2", get_rho2};
    pvd.variables += {"n.x", get_normal_x};
    pvd.variables += {"n.y", get_normal_y};

    // Создаем одномерную сетку
    Rectangle gen(0.0, 7.0, 0.0, 3.0);
    gen.set_nx(700);
    gen.set_boundaries({.left=Boundary::WALL, .right=Boundary::WALL,
                        .bottom=Boundary::WALL, .top=Boundary::WALL});

    // Создать сетку
    EuMesh mesh(gen, U);
    mesh.set_decomposition("XY");

    // Создать решатель
    MmFluid solver(mixture);
    solver.set_CFL(0.5);
    solver.set_accuracy(1);
    solver.set_method(Fluxes::CRP);
    solver.set_crp_mode(CrpMode::PLIC);
    solver.set_splitting(DirSplit::SIMPLE);

    double R_max = 1.0;
    double R_min = 0.125;
    double P_max = 1.0;
    double P_min = 0.1;

    double x_barrier = 1.0;
    double y_barrier = 1.5;

    for (auto cell: mesh) {
        Vector3d v = cell.center();

        PState z = PState::Zero();
        if (v.x() > x_barrier && v.y() < y_barrier) {
            z.density   = R_max;
            z.pressure  = P_min;
            z.mass_frac = Fractions::Pure(0);
            z.densities = ScalarSet::Pure(0, R_max);
        }
        else if (v.x() < x_barrier) {
            z.density   = R_max;
            z.pressure  = P_max;
            z.mass_frac = Fractions::Pure(1);
            z.densities = ScalarSet::Pure(1, R_max);
        }
        else {
            z.density   = R_min;
            z.pressure  = P_min;
            z.mass_frac = Fractions::Pure(2);
            z.densities = ScalarSet::Pure(2, R_min);
        }

        z.energy      = mixture.energy_rP     (z.density, z.pressure, z.mass_frac);
        z.temperature = mixture.temperature_rP(z.density, z.pressure, z.mass_frac);

        cell(U).set_state(z);
    }

    size_t n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;
    double max_time = 5.0;

    Stopwatch elapsed(true);
    while (curr_time < max_time) {
        if (curr_time >= next_write) {
            mpi::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);

            solver.interface_recovery(mesh);
            auto body_mesh0 = solver.body(mesh, 0);
            auto body_mesh1 = solver.body(mesh, 1);
            auto body_mesh2 = solver.body(mesh, 2);
            pvd_body0.save(body_mesh0, curr_time);
            pvd_body1.save(body_mesh1, curr_time);
            pvd_body2.save(body_mesh2, curr_time);

            next_write += max_time / 50;
        }

        // Точное завершение в end_time
        solver.set_max_dt(max_time - curr_time);

        // Обновляем слои
        solver.update(mesh);

        curr_time += solver.dt();
        n_step += 1;
    }
    pvd.save(mesh, max_time);

    solver.interface_recovery(mesh);
    auto body_mesh0 = solver.body(mesh, 0);
    auto body_mesh1 = solver.body(mesh, 1);
    auto body_mesh2 = solver.body(mesh, 2);
    pvd_body0.save(body_mesh0, max_time);
    pvd_body1.save(body_mesh1, max_time);
    pvd_body2.save(body_mesh2, max_time);

    elapsed.stop();

    mpi::cout << "\nElapsed:      " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    return 0;
}
