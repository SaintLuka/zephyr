/// @file Одномерные задачи газодинамики

#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/strip.h>
#include <zephyr/mesh/mesh.h>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/tests/test_1D.h>

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
double get_vfrac1(AmrStorage::Item &cell) { return cell(U).vol_frac(0); }
double get_vfrac2(AmrStorage::Item &cell) { return cell(U).vol_frac(1); }
double get_rho1(AmrStorage::Item &cell) { return cell(U).densities[0]; }
double get_rho2(AmrStorage::Item &cell) { return cell(U).densities[1]; }


int main() {
    // Тестовая задача
    //RarefiedWater test;
    //Multimat1D test(1, 1, 0);
    ToroTest test(1, true, true);

    MixturePT mixture = test.mixture_PT();

    // Файл для записи
    PvdFile pvd("test1D", "output");

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

    double curr_time = 0.0;

    pvd.variables += {"rho_exact",
                      [&test, &curr_time](const AmrStorage::Item &cell) -> double {
                          return test.density_t(cell.center, curr_time);
                      }};
    pvd.variables += {"u_exact",
                      [&test, &curr_time](const AmrStorage::Item &cell) -> double {
                          return test.velocity_t(cell.center, curr_time).x();
                      }};
    pvd.variables += {"p_exact",
                      [&test, &curr_time](const AmrStorage::Item &cell) -> double {
                          return test.pressure_t(cell.center, curr_time);
                      }};
    pvd.variables += {"e_exact",
                      [&test, &curr_time](const AmrStorage::Item &cell) -> double {
                          return test.energy_t(cell.center, curr_time);
                      }};
    pvd.variables += {"b_exact",
                      [&test, &curr_time](const AmrStorage::Item &cell) -> double {
                          return test.fractions_t(cell.center, curr_time)[0];
                      }};

    // Создаем одномерную сетку
    Strip gen(test.xmin(), test.xmax());
    gen.set_size(1000);

    // Создать сетку
    EuMesh mesh(U, &gen);

    // Создать решатель
    MmFluid solver(mixture);
    solver.set_CFL(0.5);
    solver.set_accuracy(2);
    solver.set_method(Fluxes::HLLC);

    for (auto cell: mesh) {
        Vector3d r = cell.center();
        cell(U).density  = test.density(r);
        cell(U).velocity = test.velocity(r);
        cell(U).pressure = test.pressure(r);
        cell(U).energy   = test.energy(r);

        cell(U).mass_frac = test.fractions(r);

        cell(U).densities[0] = test.fractions(r)[0] > 0.0 ? test.density(r) : NAN;
        cell(U).densities[1] = test.fractions(r)[1] > 0.0 ? test.density(r) : NAN;

        cell(U).temperature = mixture.temperature_rP(
                cell(U).density, cell(U).pressure, cell(U).mass_frac);
    }

    size_t n_step = 0;
    double next_write = 0.0;

    while (curr_time < test.max_time() && n_step < 15000000000) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);
            if (n_step < 3870000000000000000) {
                next_write += 0.0;//test.max_time() / 200;
            }
        }

        // Точное завершение в end_time
        solver.set_max_dt(test.max_time() - curr_time);

        // Обновляем слои
        solver.update(mesh);

        curr_time += solver.dt();
        n_step += 1;
    }

    // Сохранить данные как текст
    CsvFile csv("test1D.csv", 8, pvd.variables);
    csv.save(mesh);

    // Расчёт ошибок
    double r_err = 0.0, u_err = 0.0, p_err = 0.0, e_err = 0.0;
    double r_avg = 0.0, u_avg = 0.0, p_avg = 0.0, e_avg = 0.0;
    for (auto cell: mesh) {
        Vector3d r = cell.center();
        double V = cell.volume();

        r_err += V * std::abs(cell(U).density      - test.density_t(r, curr_time));
        u_err += V * std::abs(cell(U).velocity.x() - test.velocity_t(r, curr_time).x());
        p_err += V * std::abs(cell(U).pressure     - test.pressure_t(r, curr_time));
        e_err += V * std::abs(cell(U).energy       - test.energy_t(r, curr_time));

        r_avg += V * std::abs(cell(U).density     );
        u_avg += V * std::abs(cell(U).velocity.x());
        p_avg += V * std::abs(cell(U).pressure    );
        e_avg += V * std::abs(cell(U).energy      );
    }
    r_err /= r_avg;
    u_err /= u_avg;
    p_err /= p_avg;
    e_err /= e_avg;

    std::cout << "\nMean errors\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "    Density:     " << r_err << "\n";
    std::cout << "    Velocity:    " << u_err << "\n";
    std::cout << "    Pressure:    " << p_err << "\n";
    std::cout << "    Energy:      " << e_err << "\n";

    return 0;
}
