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

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::mmf;

using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;

using zephyr::utils::Stopwatch;


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


void initialize(EuMesh& mesh, const MixturePT& mixture) {
    PState z_air(
            1.0_kg_m3,          // density
            Vector3d::Zero(),   // velocity
            1.0_bar,            // pressure
            Fractions::Pure(0), // mass fractions
            mixture);

    PState z_water(
            1000.0_kg_m3,       // density
            Vector3d::Zero(),   // velocity
            1.0_bar,            // pressure
            Fractions::Pure(1), // mass fractions
            mixture);

    PState z_shock(
            1323.65_kg_m3,          // density
            {681.58_m_s, 0.0, 0.0}, // velocity
            1.9e4_bar,              // pressure
            Fractions::Pure(1),     // mass fractions
            mixture);

    double r = 3.0_mm;
    double x_bubble = 6.0_mm;
    double y_bubble = 6.0_mm;
    double x_shock  = 0.6_mm;

    Vector3d bubble_center = {x_bubble, y_bubble, 0};
    auto in_water = [bubble_center, r](const Vector3d &v) -> bool {
        return (v - bubble_center).norm() > r;
    };

    mesh.for_each([&](EuCell &cell) {
        if (cell.center().x() < x_shock) {
            cell(U).set_state(z_shock);
            return;
        }

        double vol_frac1 = cell.approx_vol_fraction(in_water);
        if (vol_frac1 == 0.0 || vol_frac1 == 1.0) {
            // Чистое вещество
            cell(U).set_state(vol_frac1 > 0.5 ? z_water : z_air);
        }
        else {
            // Граница веществ
            double vol_frac0 = cell.polygon().disk_clip_area(bubble_center, r) / cell.volume();
            double vol_frac1 = 1.0 - vol_frac0;


            mmf::PState &z0 = z_air;
            mmf::PState &z1 = z_water;

            // rho = sum a_i rho_i
            double    density   = vol_frac0 * z0.density + vol_frac1 * z1.density;
            Fractions mass_frac = {vol_frac0 * z0.density / density, vol_frac1 * z1.density / density};
            mass_frac.normalize();
            Vector3d  velocity  = mass_frac[0] * z0.velocity + mass_frac[1] * z1.velocity;
            double    pressure  = vol_frac0 * z0.pressure + vol_frac1 * z1.pressure;

            PState z(density, velocity, pressure, mass_frac, mixture);

            cell(U).set_state(z);
        }
    });
}


int main() {
    zephyr::utils::threads::off();

    // Материалы
    Eos::Ptr air = IdealGas::create("Air");
    Eos::Ptr water = StiffenedGas::create("Water");

    // Смесь
    MixturePT mixture = {air, water};

    // Файл для записи
    PvdFile pvd("BC", "output");
    PvdFile pvd_bubble("bubble", "output");

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
    Rectangle gen(0.0_cm, 1.2_cm, 0.0_cm, 1.2_cm);
    gen.set_nx(120);
    gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::WALL,
                        .bottom=Boundary::WALL, .top=Boundary::WALL});

    // Создать сетку
    EuMesh mesh(gen, U);

    // Создать решатель
    MmFluid solver(mixture);
    solver.set_CFL(0.5);
    solver.set_accuracy(1);
    solver.set_method(Fluxes::CRP);
    solver.set_crp_mode(CrpMode::PLIC);
    solver.set_splitting(DirSplit::SIMPLE);

    initialize(mesh, mixture);

    size_t n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;
    double max_time = 2.0 * 4.5_us;

    Stopwatch elapsed(true);
    while (curr_time < max_time) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << 1.0e6 * curr_time << " us\n";
            pvd.save(mesh, curr_time);

            solver.interface_recovery(mesh);
            auto bubble = solver.body(mesh, 0);
            pvd_bubble.save(bubble, curr_time);

            next_write += 0.0; //max_time / 50;
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
    auto bubble = solver.body(mesh, 0);
    pvd_bubble.save(bubble, max_time);

    elapsed.stop();

    std::cout << "\nElapsed:      " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    return 0;
}
