/// @file Одномерные задачи газодинамики.
/// Если раскомментировать #define ADAPTIVE, то задача решается
/// на двумерной адаптивной сетке.

#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/strip.h>
#include <zephyr/mesh/mesh.h>

#include <zephyr/phys/tests/test_1D.h>

#include <zephyr/math/solver/riemann.h>
#include <zephyr/math/solver/sm_fluid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/csv_file.h>

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

using zephyr::mesh::generator::Strip;
using zephyr::math::RiemannSolver;
using zephyr::mesh::EuMesh;

//#define ADAPTIVE

// Для быстрого доступа по типу
SmFluid::State U;

/// Переменные для сохранения
double get_rho(AmrStorage::Item& cell) { return cell(U).density; }
double get_u(AmrStorage::Item& cell) { return cell(U).velocity.x(); }
double get_p(AmrStorage::Item& cell) { return cell(U).pressure; }
double get_e(AmrStorage::Item& cell) { return cell(U).energy; }


int main() {
    // Тестовая задача
    SodTest test;
    //ToroTest test(1);
    //ShuOsherTest test;
    //ShockWave test(3.0, 0.1, 1.0);

    // Начальные данные
    auto init_cells = [&test](Mesh& mesh) {
        for (auto cell: mesh) {
            cell(U).density  = test.density (cell.center());
            cell(U).velocity = test.velocity(cell.center());
            cell(U).pressure = test.pressure(cell.center());
            cell(U).energy   = test.energy  (cell.center());
        }
    };

    // Уравнение состояния
    Eos::Ptr eos = test.get_eos();

    // Файл для записи
    PvdFile pvd("test1D", "output");

    double curr_time = 0.0;

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};

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
    pvd.variables += {"c",
                      [&eos](AmrStorage::Item& cell) -> double {
                          return eos->sound_speed_rP(cell(U).density, cell(U).pressure);
                      }};
    pvd.variables += {"c_exact",
                      [&eos, &test, &curr_time](const AmrStorage::Item &cell) -> double {
                          double rho = test.density_t(cell.center, curr_time);
                          double P = test.pressure_t(cell.center, curr_time);
                          return eos->sound_speed_rP(rho, P);
                      }};

    // Создаем одномерную сетку

#ifdef ADAPTIVE
    int nx = 100;
    double h = (test.xmax() - test.xmin()) / nx;

    // Костыль. Две ячейки в ширину, с одной ячейкой сейчас косяк
    // с граничными условиями
    Rectangle gen(test.xmin(), test.xmax(), -h, +h);
    gen.set_boundaries({.left   = Boundary::ZOE,  .right = Boundary::ZOE,
                        .bottom = Boundary::WALL, .top   = Boundary::WALL});
    gen.set_nx(nx);
#else
    Strip gen(test.xmin(), test.xmax());
    gen.set_boundaries({.left = Boundary::ZOE, .right = Boundary::ZOE});
    gen.set_size(500);
#endif

    // Создать сетку
    EuMesh mesh(U, &gen);

    // Создать решатель
    SmFluid solver(eos);
    solver.set_CFL(0.5);
    solver.set_accuracy(2);
    solver.set_limiter("MC");
    solver.set_method(Fluxes::HLLC);

#ifdef ADAPTIVE
    // Сеточная адаптация
    mesh.set_max_level(5);
    mesh.set_distributor(solver.distributor());

    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        init_cells(mesh);
        solver.set_flags(mesh);
        mesh.refine();
    }
#endif
    init_cells(mesh);

    size_t n_step = 0;
    double next_write = 0.0;

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

#ifdef ADAPTIVE
        solver.set_flags(mesh);
        mesh.refine();
#endif

        curr_time += solver.dt();
        n_step += 1;
    }
    pvd.save(mesh, curr_time);

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

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "\nMean errors\n";
    std::cout << "    Density:  " << r_err << "\n";
    std::cout << "    Velocity: " << u_err << "\n";
    std::cout << "    Pressure: " << p_err << "\n";
    std::cout << "    Energy:   " << e_err << "\n";

    return 0;
}
