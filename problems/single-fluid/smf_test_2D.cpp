/// @file smf_test_2D.cpp
/// @brief Двумерные газодинамические задачи, которые ставятся в прямоугольной
/// области. Используется сеточная адаптация. Можно задавать начальные условия
/// с подсеточным разрешением.

#include <iostream>
#include <iomanip>

#include <zephyr/mesh/mesh.h>

#include <zephyr/phys/tests/test_2D.h>

#include <zephyr/math/solver/sm_fluid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/csv_file.h>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/threads.h>
#include <zephyr/utils/stopwatch.h>

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

using zephyr::mesh::generator::Rectangle;
using zephyr::mesh::EuMesh;
using zephyr::math::SmFluid;
using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;


// Для быстрого доступа по типу
SmFluid::State U;

/// Переменные для сохранения
double get_lvl(AmrStorage::Item& cell) { return cell.level; }
double get_rho(AmrStorage::Item& cell) { return cell(U).density; }
double get_u(AmrStorage::Item& cell) { return cell(U).velocity.x(); }
double get_v(AmrStorage::Item& cell) { return cell(U).velocity.y(); }
double get_p(AmrStorage::Item& cell) { return cell(U).pressure; }
double get_e(AmrStorage::Item& cell) { return cell(U).energy; }


int main() {
    mpi::handler init;
    threads::on(6);

    // Тестовая задача
    Riemann2D test(6);
    //ToroTest2D test(1, 0.3 * M_PI);
    //SkewShockWave test(5.0, M_PI/6, 0.2);
    //RichtmyerMeshkov test;
    //SodTest test_1D;
    //RotatedTest test(test_1D, M_PI / 4.0);

    auto eos = test.get_eos();

    // Использовать подсеточную реконструкцию начальных данных?
    const bool simple_init = false;
    const int n = simple_init ? 0 : 20;

    // Задание начальных данных
    auto init_cells = [&test, eos, n](Mesh& mesh) {
        mesh.for_each([&](Cell& cell) {
            QState q(test.density_mean(cell, n),
                     test.momentum_mean(cell, n),
                     test.energy_mean(cell, n));
            PState z(q, *eos);
            cell(U).set_state(z);
        });
    };

    // Файл для записи
    PvdFile pvd("test2D", "output");
    pvd.unique_nodes = false;

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"v", get_v};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};

    double curr_time = 0.0;
    pvd.variables += {"rho_exact",
                      [&test, &curr_time](const AmrStorage::Item &cell) -> double {
                          return test.density_t(cell.center, curr_time);
                      }};
    pvd.variables += {"p_exact",
                      [&test, &curr_time](const AmrStorage::Item &cell) -> double {
                          return test.pressure_t(cell.center, curr_time);
                      }};

    // Генератор сетки (с граничными условиями) дает тест,
    // число ячеек можно задать
    Rectangle gen(test.xmin(), test.xmax(), test.ymin(), test.ymax());
    gen.set_boundaries(test.boundaries());
    gen.set_nx(mpi::single() ? 50 : 500);

    // Создать сетку
    EuMesh mesh(gen, U);
    mesh.set_decomposition("XY");

    // Создать решатель
    SmFluid solver(eos);
    solver.set_accuracy(2);
    solver.set_CFL(0.5);
    solver.set_limiter("MC");
    solver.set_method(Fluxes::HLLC);

    // Сеточная адаптация
    mesh.set_max_level(mpi::single() ? 3 : 0);
    mesh.set_distributor(solver.distributor());

    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        init_cells(mesh);
        solver.set_flags(mesh);
        mesh.refine();
    }
    init_cells(mesh);

    size_t n_step = 0;
    double next_write = 0.0;

    Stopwatch elapsed(true);
    while (curr_time < test.max_time()) {
        if (curr_time >= next_write) {
            mpi::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(8) << std::setprecision(3) << curr_time << "\n";

            pvd.save(mesh, curr_time);
            next_write += test.max_time() / 200;
        }
        // Точное завершение в end_time
        solver.set_max_dt(test.max_time() - curr_time);

        // Обновляем слои
        solver.update(mesh);
        solver.set_flags(mesh);
        mesh.refine();

        curr_time += solver.dt();
        n_step += 1;
    }
    pvd.save(mesh, curr_time);
    elapsed.stop();

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    // Сохранить данные как текст
    CsvFile csv("test2D.csv", 5, pvd.variables);
    csv.save(mesh);

    return 0;
}
