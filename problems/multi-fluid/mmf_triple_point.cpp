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


/// Переменные для сохранения
double get_rho(Storable<PState> part, EuCell& cell) { return cell[part].density; }
double get_vx(Storable<PState> part, EuCell& cell) { return cell[part].velocity.x(); }
double get_vy(Storable<PState> part, EuCell& cell) { return cell[part].velocity.y(); }
double get_p(Storable<PState> part, EuCell& cell) { return cell[part].pressure; }
double get_e(Storable<PState> part, EuCell& cell) { return cell[part].energy; }
double get_T(Storable<PState> part, EuCell &cell) { return cell[part].temperature; }
double get_cln(Storable<PState> part, EuCell &cell) { return cell[part].mass_frac.index(); }
double get_mfrac1(Storable<PState> part, EuCell &cell) { return cell[part].mass_frac[0]; }
double get_mfrac2(Storable<PState> part, EuCell &cell) { return cell[part].mass_frac[1]; }
double get_vfrac1(Storable<PState> part, EuCell &cell) { return cell[part].alpha(0); }
double get_vfrac2(Storable<PState> part, EuCell &cell) { return cell[part].alpha(1); }
double get_rho1(Storable<PState> part, EuCell &cell) { return cell[part].densities[0]; }
double get_rho2(Storable<PState> part, EuCell &cell) { return cell[part].densities[1]; }
double get_normal_x(Storable<VectorSet> part, EuCell &cell) { return cell[part][0].x(); }
double get_normal_y(Storable<VectorSet> part, EuCell &cell) { return cell[part][0].y(); }


int main(int argc, char** argv) {
    mpi::handler init(argc, argv);
    threads::off();

    // Материал
    Eos::Ptr sg1 = IdealGas::create(1.4, 1.0);
    Eos::Ptr sg2 = IdealGas::create(1.5, 1.0);
    Eos::Ptr sg3 = IdealGas::create(1.5, 1.0);

    // Формальная смесь
    MixturePT mixture = {sg1, sg2, sg3};

    // Создаем одномерную сетку
    Rectangle gen(0.0, 7.0, 0.0, 3.0);
    gen.set_nx(700);
    gen.set_boundaries({.left=Boundary::WALL, .right=Boundary::WALL,
                        .bottom=Boundary::WALL, .top=Boundary::WALL});

    // Создать сетку
    EuMesh mesh(gen);
    mesh.set_decomposition("XY");

    // Создать решатель
    MmFluid solver(mixture);
    solver.set_CFL(0.5);
    solver.set_accuracy(1);
    solver.set_method(Fluxes::CRP);
    solver.set_crp_mode(CrpMode::PLIC);
    solver.set_splitting(DirSplit::SIMPLE);
    solver.add_types(mesh);

    // Файл для записи
    PvdFile pvd("TP", "output");
    PvdFile pvd_body0("body0", "output");
    PvdFile pvd_body1("body1", "output");
    PvdFile pvd_body2("body2", "output");

    // Переменные для сохранения
    pvd.variables += {"cln", std::bind(get_cln, solver.part.init, std::placeholders::_1) };
    pvd.variables += {"rho", std::bind(get_rho, solver.part.init, std::placeholders::_1)};
    pvd.variables += {"vx", std::bind(get_vx, solver.part.init, std::placeholders::_1)};
    pvd.variables += {"vy", std::bind(get_vy, solver.part.init, std::placeholders::_1)};
    pvd.variables += {"e",  std::bind(get_e, solver.part.init, std::placeholders::_1)};
    pvd.variables += {"P",  std::bind(get_p, solver.part.init, std::placeholders::_1)};
    pvd.variables += {"T",  std::bind(get_T, solver.part.init, std::placeholders::_1)};
    pvd.variables += {"b1", std::bind(get_mfrac1, solver.part.init, std::placeholders::_1)};
    pvd.variables += {"b2", std::bind(get_mfrac2, solver.part.init, std::placeholders::_1)};
    pvd.variables += {"a1", std::bind(get_vfrac1, solver.part.init, std::placeholders::_1)};
    pvd.variables += {"a2", std::bind(get_vfrac2, solver.part.init, std::placeholders::_1)};
    pvd.variables += {"rho1", std::bind(get_rho1, solver.part.init, std::placeholders::_1)};
    pvd.variables += {"rho2", std::bind(get_rho2, solver.part.init, std::placeholders::_1)};
    pvd.variables += {"n.x", std::bind(get_normal_x, solver.part.n, std::placeholders::_1)};
    pvd.variables += {"n.y", std::bind(get_normal_y, solver.part.n, std::placeholders::_1)};

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
