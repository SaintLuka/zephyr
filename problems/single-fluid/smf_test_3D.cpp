/// @file smf_test_3D.cpp
/// @brief Сферический взрыв на трехмерной декартовой сетке.


#include <iostream>
#include <iomanip>

#include <zephyr/phys/tests/test_3D.h>

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

using zephyr::mesh::generator::Cuboid;
using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::math::SmFluid;
using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;

// Критерий адаптации подобран под задачу.
// Адаптация ячеек с плотностью выше 1.5.
void set_flags(EuMesh &mesh, Storable<PState> z) {
    if (!mesh.adaptive()) return;

    mesh.for_each([z](EuCell cell) {
        const double threshold = 1.5;

        if (cell(z).density > threshold) {
            cell.set_flag(1);
            return;
        }
        for (auto face: cell.faces()) {
            if (face.neib(z).density > threshold) {
                cell.set_flag(1);
                return;
            }
        }
        cell.set_flag(-1);
    });
}

int main(int argc, char** argv) {
    mpi::handler handler(argc, argv);
    threads::init(argc, argv);
    threads::info();

    // Тестовая задача
    SedovBlast3D test({.gamma=1.4, .rho0=1.0, .E=1.0});

    // Уравнение состояния
    Eos::Ptr eos = test.get_eos();

    // Генератор сетки
    Cuboid gen = Cuboid(test.xmin(), test.xmax(),
                        test.ymin(), test.ymax(),
                        test.zmin(), test.zmax());
    gen.set_nx(10);
    gen.set_boundaries(test.boundaries());

    // Создать сетку
    EuMesh mesh(gen);

    // Создать решатель
    SmFluid solver(eos);
    solver.set_accuracy(1);
    solver.set_CFL(0.5);
    solver.set_limiter("minmod");
    solver.set_method(Fluxes::HLL);

    // Добавляем типы на сетку, выбираем основной слой
    auto data = solver.add_types(mesh);
    auto z = data.init;

    // Настройка сетки
    //mesh.set_decomposition("XYZ");
    //mesh.block_sorting(4, 4, 8);
    mesh.set_max_level(4);
    mesh.set_distributor(solver.distributor());

    // Начальные данные
    auto init_cells = [&](EuMesh& mesh) {
        mesh.for_each([&](EuCell& cell) {
            Vector3d r = cell.center();
            cell(z).density  = test.density (r);
            cell(z).velocity = test.velocity(r);
            cell(z).pressure = std::max(test.pressure(r), 1.0e-3);
            cell(z).energy   = eos->energy_rP(cell(z).density, cell(z).pressure);
        });
    };

    // Файл для записи
    PvdFile pvd("Sedov", "output");

    size_t n_step = 0;
    double curr_time = test.init_time;
    double next_write = test.init_time;

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables += {"rho", [z](EuCell& cell) -> double { return cell(z).density; }};
    pvd.variables += {"vr",  [z](EuCell& cell) -> double { return cell(z).velocity.norm(); }};
    pvd.variables += {"p",   [z](EuCell& cell) -> double { return cell(z).pressure; }};
    pvd.variables += {"e",   [z](EuCell& cell) -> double { return cell(z).energy; }};
    pvd.variables += {"rho_exact",
                      [&test, &curr_time](const EuCell &cell) -> double {
                          return test.density_t(cell.center(), curr_time);
                      }};
    pvd.variables += {"vr_exact",
                      [&test, &curr_time](const EuCell  &cell) -> double {
                          return test.velocity_t(cell.center(), curr_time).norm();
                      }};
    pvd.variables += {"p_exact",
                      [&test, &curr_time](const EuCell &cell) -> double {
                          return test.pressure_t(cell.center(), curr_time);
                      }};
    pvd.variables += {"e_exact",
                      [&test, &curr_time](const EuCell &cell) -> double {
                          return test.energy_t(cell.center(), curr_time);
                      }};

    // Инициализация начальными данными
    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        init_cells(mesh);
        set_flags(mesh, z);
        mesh.refine();
    }
    init_cells(mesh);

    Stopwatch sw_update;
    Stopwatch sw_flags;
    Stopwatch sw_grad;
    Stopwatch sw_refine;

    Stopwatch elapsed(true);
    while (curr_time < test.max_time()) {
        if (curr_time >= next_write) {
            mpi::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(10) << std::setprecision(5) << curr_time << "\n";

            pvd.save(mesh, curr_time);
            next_write += test.max_time() / 200;
        }

        // Точное завершение в end_time
        solver.set_max_dt(test.max_time() - curr_time);

        // Обновляем слои
        sw_update.resume();
        solver.update(mesh);
        sw_update.stop();

        sw_flags.resume();
        set_flags(mesh, z);
        sw_flags.stop();

        // Для переноса по градиентам
        sw_grad.resume();
        solver.compute_grad(mesh);
        sw_grad.stop();

        sw_refine.resume();
        mesh.refine();
        sw_refine.stop();

        curr_time += solver.dt();
        n_step += 1;
    }
    elapsed.stop();
    pvd.save(mesh, 1.0);

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    mpi::cout << "  Update time:  " << sw_update.extended_time()
              << " ( " << sw_update.milliseconds() << " ms)\n";
    mpi::cout << "  Flags  time:  " << sw_flags.extended_time()
              << " ( " << sw_flags.milliseconds() << " ms)\n";
    mpi::cout << "  Gradient   :  " << sw_grad.extended_time()
              << " ( " << sw_grad.milliseconds() << " ms)\n";
    mpi::cout << "  Refine time:  " << sw_refine.extended_time()
              << " ( " << sw_refine.milliseconds() << " ms)\n";

    return 0;
}