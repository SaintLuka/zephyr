/// @file fb_test_1D.cpp
/// @brief Одномерные задачи газодинамики. Если раскомментировать #define ADAPTIVE,
/// то задача решается на двумерной адаптивной сетке.

#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/strip.h>

#include <zephyr/phys/tests/test_1D.h>

#include <zephyr/math/solver/riemann.h>
#include <zephyr/math/solver/free_boundary.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/csv_file.h>

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

using zephyr::geom::generator::Strip;
using zephyr::math::RiemannSolver;
using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::math::FreeBoundary;
using zephyr::utils::threads;
using zephyr::utils::mpi;

#define ADAPTIVE

int main() {
    mpi::handler init;
    threads::on();

    // Тестовая задача
    //SodTest test;
    //ToroTest test(1);
    ShuOsherTest test;
    //ShockWave test(3.0, 0.1, 1.0);

    // Уравнение состояния
    Eos::Ptr eos = test.get_eos();

    // Создаем одномерную сетку
#ifdef ADAPTIVE
    int nx = 100;
    double h = 0.01 * (test.xmax() - test.xmin());

    Rectangle gen(test.xmin(), test.xmax(), -h, +h);
    gen.set_boundaries({.left   = Boundary::ZOE,  .right = Boundary::ZOE,
                        .bottom = Boundary::WALL, .top   = Boundary::WALL});
    gen.set_sizes(nx, 1);
#else
    Strip gen(test.xmin(), test.xmax());
    gen.set_boundaries({.left = Boundary::ZOE, .right = Boundary::ZOE});
    gen.set_size(500);
#endif

    // Создать сетку
    EuMesh mesh(gen);

    // Создать решатель
    FreeBoundary solver(eos);
    solver.set_CFL(0.5);
    solver.set_accuracy(2);
    solver.set_limiter("MC");
    solver.set_method(Fluxes::HLLC_M);

    // Добавляем типы на сетку, выбираем основной слой
    auto data = solver.add_types(mesh);
    auto z = data.init;

    // Файл для записи
    PvdFile pvd("test1D", "D:\\main_project");

    double curr_time = 0.0;

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables += {"density",  [z](EuCell& cell) -> double { return cell(z).density; }};
    pvd.variables += {"vel.x",    [z](EuCell& cell) -> double { return cell(z).velocity.x(); }};
    pvd.variables += {"vel.y",    [z](EuCell& cell) -> double { return cell(z).velocity.y(); }};
    pvd.variables += {"pressure", [z](EuCell& cell) -> double { return cell(z).pressure; }};
    pvd.variables += {"energy",   [z](EuCell& cell) -> double { return cell(z).energy; }};

    pvd.variables += {"exact.dens",
                      [&test, &curr_time](const EuCell &cell) -> double {
                          return test.density_t(cell.center(), curr_time);
                      }};
    pvd.variables += {"exact.velocity",
                      [&test, &curr_time](const EuCell &cell) -> double {
                          return test.velocity_t(cell.center(), curr_time).x();
                      }};
    pvd.variables += {"exact.pres",
                      [&test, &curr_time](const EuCell &cell) -> double {
                          return test.pressure_t(cell.center(), curr_time);
                      }};
    pvd.variables += {"exact.energy",
                      [&test, &curr_time](const EuCell &cell) -> double {
                          return test.energy_t(cell.center(), curr_time);
                      }};
    pvd.variables += {"sound_speed",
                      [&eos, z](const EuCell & cell) -> double {
                          return eos->sound_speed_rP(cell(z).density, cell(z).pressure);
                      }};
    pvd.variables += {"exact.sound",
                      [&eos, &test, &curr_time](const EuCell &cell) -> double {
                          double rho = test.density_t(cell.center(), curr_time);
                          double P = test.pressure_t(cell.center(), curr_time);
                          return eos->sound_speed_rP(rho, P);
                      }};

    // Начальные данные
    auto init_cells = [&test, z](EuMesh& mesh) {
        for (auto cell: mesh) {
            cell(z).density  = test.density (cell.center());
            cell(z).velocity = test.velocity(cell.center());
            cell(z).pressure = test.pressure(cell.center());
            cell(z).energy   = test.energy  (cell.center());
        }
    };

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
    csv.save(mesh.locals());

    // Расчёт ошибок
    double r_err = 0.0, u_err = 0.0, p_err = 0.0, e_err = 0.0;
    double r_avg = 0.0, u_avg = 0.0, p_avg = 0.0, e_avg = 0.0;
    for (auto cell: mesh) {
        Vector3d r = cell.center();
        double V = cell.volume();

        r_err += V * std::abs(cell(z).density      - test.density_t(r, curr_time));
        u_err += V * std::abs(cell(z).velocity.x() - test.velocity_t(r, curr_time).x());
        p_err += V * std::abs(cell(z).pressure     - test.pressure_t(r, curr_time));
        e_err += V * std::abs(cell(z).energy       - test.energy_t(r, curr_time));

        r_avg += V * std::abs(cell(z).density     );
        u_avg += V * std::abs(cell(z).velocity.x());
        p_avg += V * std::abs(cell(z).pressure    );
        e_avg += V * std::abs(cell(z).energy      );
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
