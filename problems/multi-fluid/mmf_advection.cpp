/// @file Двумерная задача переноса

#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/strip.h>
#include <zephyr/mesh/mesh.h>

#include <zephyr/phys/tests/toro.h>
#include <zephyr/phys/tests/rarefied_water.h>
#include <zephyr/phys/tests/multimat_1D.h>

#include <zephyr/math/solver/riemann.h>
#include <zephyr/math/solver/mm_fluid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/csv_file.h>

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::mmf;

using zephyr::mesh::generator::Strip;
using zephyr::math::RiemannSolver;
using zephyr::mesh::EuMesh;


// Для быстрого доступа по типу
MmFluid::State U;

/// Переменные для сохранения
double get_rho(AmrStorage::Item& cell) { return cell(U).density; }
double get_u(AmrStorage::Item& cell) { return cell(U).velocity.x(); }
double get_p(AmrStorage::Item& cell) { return cell(U).pressure; }
double get_e(AmrStorage::Item& cell) { return cell(U).energy; }
double get_T(AmrStorage::Item &cell) { return cell(U).temperature; }
double get_cln(AmrStorage::Item &cell) { return cell(U).mass_frac.index(); }
double get_mfrac1(AmrStorage::Item &cell) { return cell(U).mass_frac[0]; }
double get_mfrac2(AmrStorage::Item &cell) { return cell(U).mass_frac[1]; }
double get_vfrac1(AmrStorage::Item &cell) { return cell(U).vol_frac[0]; }
double get_vfrac2(AmrStorage::Item &cell) { return cell(U).vol_frac[1]; }
double get_rho1(AmrStorage::Item &cell) { return cell(U).get_state().true_density(0); }
double get_rho2(AmrStorage::Item &cell) { return cell(U).get_state().true_density(1); }
double get_normal_x(AmrStorage::Item &cell) { return cell(U).n[0].x(); }
double get_normal_y(AmrStorage::Item &cell) { return cell(U).n[0].y(); }


int main() {
    // Материал
    double gamma = 1.4;
    Eos::Ptr sg1 = IdealGas::create(gamma, 1.0);
    Eos::Ptr sg2 = IdealGas::create(gamma, 1.0);

    // Формальная смесь
    MixturePT mixture = {sg1, sg2};

    // Файл для записи
    PvdFile pvd("Transfer2D", "output");

    // Переменные для сохранения
    pvd.variables += {"cln", get_cln};
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
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
    Rectangle gen(0.0, 1.0, 0.0, 0.7);
    gen.set_nx(200);
    gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
                        .bottom=Boundary::ZOE, .top=Boundary::ZOE});

    // Создать сетку
    EuMesh mesh(U, &gen);

    // Создать решатель
    MmFluid solver(mixture);
    solver.set_CFL(0.5);
    solver.set_accuracy(1);
    solver.set_method(Fluxes::CRP);
    solver.set_crp_mode(CrpMode::PLIC);
    solver.set_splitting(DirSplit::SIMPLE);

    for (auto cell: mesh) {
        Vector3d r = cell.center();

        bool in = std::abs(r.x() - 0.15) < 0.1 &&
                  std::abs(r.y() - 0.50) < 0.1;

        cell(U).velocity    = {0.7, -0.35, 0.0};
        cell(U).density     = 1.0;
        cell(U).pressure    = 1.0 / gamma;
        //cell(U).energy      = 1.0 / (gamma * (gamma - 1.0));
        //cell(U).temperature = 1.0;

        /*
        if (r.x() < 0.1) {
            cell(U).mass_frac[0]  = 0.0;
        }
        if (std::abs(r.x() - 0.2) < 0.1) {
            cell(U).mass_frac[0]  = (r.x() - 0.1) / 0.2;
        }
        if (r.x() > 0.3) {
            cell(U).mass_frac[0]  = 1.0;
        }
         */
        cell(U).mass_frac[0] = in ? 1.0 : 0.0;
        cell(U).mass_frac[1] = 1.0 - cell(U).mass_frac[0] ;

        cell(U).vol_frac[0] = cell(U).mass_frac[0];
        cell(U).vol_frac[1] = cell(U).mass_frac[1];

        /*
        cell(U).density = 1.0 / mixture.volume_PT(
                cell(U).pressure, cell(U).temperature, cell(U).mass_frac);
        cell(U).energy = mixture.energy_PT(
                cell(U).pressure, cell(U).temperature, cell(U).mass_frac);
                */
        cell(U).energy      = mixture.energy_rP(
                cell(U).density, cell(U).pressure, cell(U).mass_frac);
        cell(U).temperature = mixture.temperature_rP(
                cell(U).density, cell(U).pressure, cell(U).mass_frac);

    }

    size_t n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;
    double max_time = 1.0;

    while (curr_time < max_time) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);
            next_write += max_time / 200;
        }

        // Точное завершение в end_time
        solver.set_max_dt(max_time - curr_time);

        // Обновляем слои
        solver.update(mesh);

        curr_time += solver.dt();
        n_step += 1;
    }
    pvd.save(mesh, max_time);

    return 0;
}
