/// @file smf_test_2D.cpp

#include <iostream>
#include <iomanip>

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

using zephyr::geom::generator::Rectangle;
using zephyr::mesh::Mesh;
using zephyr::mesh::Cell;
using zephyr::math::SmFluid;
using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;


int main(int argc, char** argv) {
    mpi::handler handler(argc, argv);
    threads::init(argc, argv);
    threads::info();
    threads::on(7);

    // Файл для записи
    PvdFile pvd("mesh", "output");
    pvd.options.unique_nodes = true;

    Stopwatch elapsed(true);
    // case 5, 13, 14, 15 не захватывает ни в какую небольшие волны
    for (int test_case = 1; test_case <= 19; ++test_case) {
        // Тестовая задача
        Riemann2D test(test_case);
        auto eos = test.get_eos();

        // Генератор сетки (с граничными условиями) дает тест
        Rectangle gen(test.xmin(), test.xmax(), test.ymin(), test.ymax());
        gen.set_boundaries(test.boundaries());
        gen.set_nx(200);

        // Создать сетку
        Mesh mesh(gen);

        // Создать и настроить решатель
        SmFluid solver(eos);
        solver.set_accuracy(2);
        solver.set_CFL(0.5);
        solver.set_limiter(test_case == 11 ? "van Albada" : "MC");
        solver.set_method(Fluxes::HLLC);
        solver.set_criterion(SlopeCriterion::DensPres(0.04, 0.04));

        // Добавляем типы на сетку, выбираем основной слой
        auto data = solver.add_types(mesh);
        auto z = data.init;

        // Настройка сетки
        mesh.set_decomposition("XY");
        mesh.set_max_level(3);
        mesh.set_distributor(solver.distributor());

        // Задание начальных данных
        auto init_cells = [&test, z, eos](Mesh& mesh) {
            mesh.for_each([&](Cell& cell) {
                Vector3d r = cell.center();
                cell[z].density  = test.density(r);
                cell[z].velocity = test.velocity(r);
                cell[z].pressure = test.pressure(r);
                cell[z].energy   = eos->energy_rP(cell[z].density, cell[z].pressure);
            });
        };

        // Переменные для сохранения
        pvd.variables = {"level"};
        pvd.variables += {"density",  [z](Cell& cell) -> double { return cell[z].density; }};
        pvd.variables += {"vel.x",    [z](Cell& cell) -> double { return cell[z].velocity.x(); }};
        pvd.variables += {"vel.y",    [z](Cell& cell) -> double { return cell[z].velocity.y(); }};
        pvd.variables += {"pressure", [z](Cell& cell) -> double { return cell[z].pressure; }};
        pvd.variables += {"energy",   [z](Cell& cell) -> double { return cell[z].energy; }};

        double curr_time = 0.0;

        // Инициализация начальными данными
        for (int k = 0; k < mesh.max_level() + 3; ++k) {
            init_cells(mesh);
            solver.set_flags(mesh);
            mesh.refine();
        }
        init_cells(mesh);

        size_t n_step = 0;
        double next_write = 0.0;

        mpi::cout << "Run test #" << test_case << "\n";
        while (curr_time < test.max_time()) {
            // Load balancing
            if (n_step % 20 == 0) {
                mesh.balancing(mesh.n_cells());
                mesh.redistribute(z);
            }

            // Print
            if (curr_time >= next_write) {
                mpi::cout << "\tStep: " << std::setw(6) << n_step << ";"
                          << "\tTime: " << std::setw(8) << std::setprecision(3) << curr_time << "\n";
                next_write += test.max_time() / 20;
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
        pvd.save(mesh, test_case);
    }
    elapsed.stop();

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";
    return 0;
}
