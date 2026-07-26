/// @file swe_lake.cpp
/// @brief Статичный двумерный тест озеро в покое.
/// Well-balanced схема демонстрирует отсутствие осцилляций.

#include <iomanip>

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

using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;

int main(int argc, char** argv) {
    mpi::handler handler(argc, argv);
    threads::init(argc, argv);
    threads::info();

	// Генератор сетки
	generator::Rectangle gen(-0.5, 0.5, -0.5, 0.5);
	gen.set_boundaries({.left = Boundary::WALL, .right = Boundary::WALL,
						.bottom = Boundary::WALL, .top=Boundary::WALL});
	gen.set_nx(200);

    // Создать сетку
    EuMesh mesh(gen);

    // Создать и настроить решатель
    SwSolver solver(PitBed::create(0.4, -0.4, 0.0));
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
	pvd.unique_nodes = true;

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables += {"surf",  [z,zb](EuCell& cell) -> double { return cell[z].surf(cell[zb]); }};
    pvd.variables += {"bed",   [zb](EuCell& cell) -> double { return cell[zb]; }};
    pvd.variables += {"depth", [z](EuCell& cell) -> double { return cell[z].h(); }};
    pvd.variables += {"vel.x", [z](EuCell& cell) -> double { return cell[z].u(); }};
    pvd.variables += {"vel.y", [z](EuCell& cell) -> double { return cell[z].v(); }};

	// Задание начальных данных
	auto init_cells = [&]() {
		mesh.for_each([&](EuCell& cell) {
			double level = -0.1;
		    cell[z].depth = std::max(0.0, level - cell[zb]);
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
    double curr_time = 0.0;
    double max_time = 5.0;

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
    return 0;
}
