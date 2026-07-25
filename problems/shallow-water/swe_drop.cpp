/// @file swe_drop.cpp
/// @brief ???

#include <iomanip>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/tests/swe/dam_break.h>
#include <zephyr/phys/matter/eos/mie_gruneisen.h>
#include <zephyr/geom/generator/strip.h>
#include <zephyr/geom/generator/rectangle.h>

#include <zephyr/math/solver/sw_solver.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/csv_file.h>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/threads.h>
#include <zephyr/utils/stopwatch.h>

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::math::SwSolver;
using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;

int main(int argc, char** argv) {
    mpi::handler handler(argc, argv);
    threads::init(argc, argv);
    threads::info();
    threads::off();

    // Генератор сетки
	// generator::Strip gen(-0.5, 0.5);
	// gen.set_boundaries({.left = Boundary::ZOE, .right = Boundary::ZOE});
	// gen.set_nx(10000);

	generator::Rectangle gen(-0.5, 0.5, -0.5, 0.5);
	gen.set_boundaries({.left = Boundary::WALL, .right = Boundary::WALL,
						.bottom = Boundary::WALL, .top=Boundary::WALL});
	gen.set_nx(100);

    // Создать сетку
    EuMesh mesh(gen);

    MieGruneisen::Ptr eos = MieGruneisen::create("Fe");

    // Создать и настроить решатель
    SwSolver solver(-0.3);
    solver.set_accuracy(1);
    solver.set_CFL(0.5);
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
    DamBreak test(2.0, 1.0);

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables += {"eta", [z](EuCell& cell) -> double { return cell[z].level; }};
    pvd.variables += {"bed", [zb](EuCell& cell) -> double { return cell[zb]; }};
    pvd.variables += {"h",   [z, zb](EuCell& cell) -> double { return cell[z].depth(cell[zb]); }};
    pvd.variables += {"u",   [z](EuCell& cell) -> double { return cell[z].velocity.x(); }};
    pvd.variables += {"v",   [z](EuCell& cell) -> double { return cell[z].velocity.y(); }};

    pvd.variables += {"eta.exact", [test, &curr_time](EuCell& cell) -> double { return test.level(cell.x(), curr_time); }};
    pvd.variables += {"h.exact",   [test, &curr_time](EuCell& cell) -> double { return test.depth(cell.x(), curr_time); }};
    pvd.variables += {"u.exact",   [test, &curr_time](EuCell& cell) -> double { return test.speed(cell.x(), curr_time); }};

	// Задание начальных данных
	auto init_cells = [&]() {
		mesh.for_each([&](EuCell& cell) {
		    cell[z].level = cell.x() < 0.0 ? 0.0 : 1.0; //0.2 * std::exp(-100.0 * cell.center().squaredNorm());
		    cell[z].velocity = Vector2d::Zero();
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
    double max_time = 0.05;

    Stopwatch elapsed(true);
    while (curr_time < max_time) {
        if (curr_time >= next_write) {
            mpi::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(8) << std::setprecision(3) << curr_time << "\n";

            pvd.save(mesh, curr_time);
            next_write += max_time / 50;
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
    elapsed.stop();

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    return 0;
}
