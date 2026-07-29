/// @file swe_test_1D.cpp
/// @brief Одномерные тесты, часто с аналитическим решением.

#include <iomanip>

#include <zephyr/geom/generator/strip.h>
#include <zephyr/phys/tests/swe/dam_break.h>
#include <zephyr/phys/tests/swe/thacker.h>
#include <zephyr/phys/tests/swe/step.h>

#include <zephyr/math/solver/sw_solver.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/csv_file.h>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/threads.h>
#include <zephyr/utils/stopwatch.h>

using namespace zephyr::io;
using namespace zephyr::phys::swe;
using namespace zephyr::math;

using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::math::SwSolver;
using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;

int main() {
    threads::off();

    //DamBreak test(2);
	//Thacker1D test;
	Step test;

	// Генератор сетки
	generator::Strip gen(test.x_min(), test.x_max());
	gen.set_boundaries({.left = Boundary::ZOE, .right = Boundary::ZOE});
	gen.set_nx(500);

    // Создать сетку
    EuMesh mesh(gen);

    // Создать и настроить решатель
    SwSolver solver(test.topography());
    solver.set_accuracy(1);
    solver.set_CFL(0.2);
    solver.set_limiter("MC");
    solver.set_method(Fluxes::HLL);

    // Добавляем типы на сетку, выбираем основной слой
    auto data = solver.add_types(mesh);
    auto z = data.init;
	auto zb = data.bed;

    // Настройка сетки
    mesh.set_decomposition("XY");
    mesh.set_max_level(0);
    mesh.set_distributor(solver.distributor());

    // Файл для записи
    PvdFile pvd("mesh", "output");

    double curr_time = 0.0;

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables += {"surf",  [z,zb](EuCell& cell) -> double { return cell[z].surf(cell[zb]); }};
    pvd.variables += {"bed",   [zb](EuCell& cell) -> double { return cell[zb]; }};
    pvd.variables += {"depth", [z](EuCell& cell) -> double { return cell[z].h(); }};
    pvd.variables += {"vel.x", [z](EuCell& cell) -> double { return cell[z].u(); }};
    pvd.variables += {"vel.y", [z](EuCell& cell) -> double { return cell[z].v(); }};

    pvd.variables += {"exact.surf",  [test, &curr_time](EuCell& cell) -> double { return test.level(cell.x(), curr_time); }};
    pvd.variables += {"exact.depth", [test, &curr_time](EuCell& cell) -> double { return test.depth(cell.x(), curr_time); }};
    pvd.variables += {"exact.vel",   [test, &curr_time](EuCell& cell) -> double { return test.speed(cell.x(), curr_time); }};

	// Задание начальных данных
	auto init_cells = [&]() {
		mesh.for_each([&](EuCell& cell) {
			cell[zb] = test.bed(cell.x());
			cell[z].depth = test.depth(cell.x(), curr_time);
			cell[z].velocity.x() = test.speed(cell.x(), curr_time);
			cell[z].velocity.y() = 0.0;
	    });
	};

    // Инициализация начальными данными
    for (int k = 0; mesh.adaptive() && (k < mesh.max_level() + 3); ++k) {
        init_cells();
        solver.set_flags(mesh);
        mesh.refine();
    }
    init_cells();

    size_t n_step = 0;
    double next_write = 0.0;
    double max_time = test.max_time();

    Stopwatch elapsed(true);
    while (curr_time < max_time) {
        if (curr_time >= next_write) {
            mpi::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(8) << std::setprecision(3) << curr_time << "\n";

            pvd.save(mesh, curr_time);
            next_write += max_time / 100;
        }

        // Точное завершение в end_time
        solver.set_max_dt(max_time - curr_time);

        // Обновляем слои
        solver.update(mesh);
        solver.set_flags(mesh);
        mesh.refine();

        curr_time += solver.dt();
        n_step += 1;
    }
    pvd.save(mesh, curr_time);
    return 0;
}
