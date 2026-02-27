/// @file smf_skew_wave.cpp
/// @brief ???

#include <iomanip>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/matter/eos/mie_gruneisen.h>
#include <zephyr/geom/generator/rectangle.h>

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

using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::math::SmFluid;
using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;

// Усредненный консервативный вектор состояния в ячейке
QState mean(EuCell& cell, const std::function<QState(const Vector3d&)>& get_state, int n) {
	auto density    = [&get_state](const Vector3d& r) -> double { return get_state(r).density; };
	auto momentum_x = [&get_state](const Vector3d& r) -> double { return get_state(r).momentum.x(); };
	auto momentum_y = [&get_state](const Vector3d& r) -> double { return get_state(r).momentum.y(); };
	auto momentum_z = [&get_state](const Vector3d& r) -> double { return get_state(r).momentum.z(); };
	auto energy     = [&get_state](const Vector3d& r) -> double { return get_state(r).energy; };

	if (n < 2) {
		return get_state(cell.center());
	}
	if (cell.const_function(density) &&
		cell.const_function(momentum_x) &&
		cell.const_function(momentum_y) &&
		cell.const_function(momentum_z) &&
		cell.const_function(energy)) {
		return get_state(cell.center());
	}

	double V = cell.volume();

	QState q;
	q.density = cell.integrate_low(density, n) / V;
	q.momentum.x() = cell.integrate_low(momentum_x, n) / V;
	q.momentum.y() = cell.integrate_low(momentum_y, n) / V;
	q.momentum.z() = cell.integrate_low(momentum_z, n) / V;

	return q;
}

int main(int argc, char** argv) {
    mpi::handler handler(argc, argv);
    threads::init(argc, argv);
    threads::info();

    // Генератор сетки
	double Ly = 0.1_cm;
	double v2x = 3.145830462176495_km; // Скорость в ЛабСО из MatLab
	double v1x = 5.0_km; //Х-скорость налёта среды на УВ
	double vy = 6.0_km; //У-скорость в ЛабСО, изменяемый параметр
	double v2 = sqrt(pow(v2x,2)+pow(vy,2)); // Модуль скорсоти перед фронтом УВ в ЛабСО
	double alpha1 = atan(vy/v1x); //Угол между v1 и осью Х в ЛабСО
	double alpha2 = atan(vy/v2x); //Угол между v2 и осью Х в ЛабСО
	double nx = -cos(alpha1); // Х-компонента нормали фронта УВ в СО, сонаправленной с v1
	double ny = sin(alpha1);// У-компонента нормали фронта УВ в СО, сонаправленной с v1
	double v2x_newSO = v2*cos(alpha2-alpha1); // Х-компонента скорости за фронтом УВ в СО, сонаправленной с v1
	double v2y_newSO = v2*sin(alpha2-alpha1); // У-компонента скорости за фронтом УВ в СО, сонаправленной с v1
	double ky = 2.0*M_PI / Ly;
	double kx = 3.0*ky;

    generator::Rectangle gen(-0.5_cm, 0.5_cm, -5*Ly, 5*Ly);
    gen.set_boundaries({.left = Boundary::ZOE, .right = Boundary::ZOE,
                        .bottom = Boundary::ZOE, .top=Boundary::ZOE});
    gen.set_nx(1000);

    // Создать сетку
    EuMesh mesh(gen);

    MieGruneisen::Ptr eos = MieGruneisen::create("Fe");

    // Создать и настроить решатель
    SmFluid solver(eos);
    solver.set_accuracy(2);
    solver.set_CFL(0.5);
    solver.set_limiter("MC");
    solver.set_method(Fluxes::HLLC);

    // Добавляем типы на сетку, выбираем основной слой
    auto data = solver.add_types(mesh);
    auto z = data.init;

    // Настройка сетки
    mesh.set_decomposition("XY");
    mesh.set_max_level(0);
    mesh.set_distributor(solver.distributor());

    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables += {"density",  [z](EuCell& cell) -> double { return cell[z].density; }};
    pvd.variables += {"vel.x",    [z](EuCell& cell) -> double { return cell[z].velocity.x(); }};
    pvd.variables += {"vel.y",    [z](EuCell& cell) -> double { return cell[z].velocity.y(); }};
    pvd.variables += {"pressure", [z](EuCell& cell) -> double { return cell[z].pressure; }};
    pvd.variables += {"energy",   [z](EuCell& cell) -> double { return cell[z].energy; }};

	// Выдает консервативный вектор состояния
	auto get_state = [&](const Vector3d& r) -> QState {
		PState st;
		if (r.x()*nx + r.y()*ny < 0.0) {
			// За фронтом УВ
			st.density  = 9.995414342253367_g_cm3;
			st.velocity.x() = v2x_newSO-sqrt(pow(v1x,2)+pow(vy,2));
			st.velocity.y() = v2y_newSO;
			st.pressure = 83.028520532665250_GPa;
		}
		else {
			// Перед фронтом УВ
			st.density  = eos->density();
			st.velocity.x() = 0;
			st.velocity.y() = 0;
			st.pressure = 1.0_bar;

			double rho0 = eos->density();

			st.density = (1.0 + 0.05 * std::cos(kx * r.x() + ky * (r.y() - 0.25 * Ly))) * rho0*0.9;
		}
		st.energy = eos->energy_rP(st.density, st.pressure);
		return QState(st);
	};

	// Задание начальных данных
	auto init_cells=[&]() {
		mesh.for_each([&](EuCell& cell) {
			QState q = mean(cell, get_state, 5);
			cell[z] = PState(q, *eos);
	    });
	};


    // Инициализация начальными данными
    for (int k = 0; mesh.adaptive() && k < mesh.max_level() + 3; ++k) {
        init_cells();
        solver.set_flags(mesh);
        mesh.refine();
    }
    init_cells();

    size_t n_step = 0;
    double next_write = 0.0;
    double curr_time = 0.0;
    double max_time = 1.5_us;

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
