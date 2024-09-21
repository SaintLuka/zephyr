/// @file Двумерные газодинамические задачи, которые ставятся в прямоугольной области.
/// Используется сеточная адаптация. Можно задавать начальные условия с подсеточным разрешением.

#include <iostream>
#include <iomanip>

#include <zephyr/mesh/mesh.h>

#include <zephyr/phys/tests/toro.h>
#include <zephyr/phys/tests/sedov.h>
#include <zephyr/phys/tests/test_2D.h>
#include <zephyr/phys/tests/shu_osher.h>
#include <zephyr/phys/tests/shock_wave.h>
#include <zephyr/phys/tests/richtmyer_meshkov.h>

#include <zephyr/math/solver/sm_fluid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/csv_file.h>

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

using zephyr::mesh::generator::Rectangle;
using zephyr::mesh::EuMesh;
using zephyr::math::SmFluid;
using zephyr::utils::threads;


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
    threads::on(6);

    // Тестовая задача
    //Test2D test(6);
    SedovBlast test;
    //ShuOsherTest test;
    //ToroTest2D test(1, 0.3 * M_PI);
    //SkewShockWave test(5.0, M_PI/6, 0.2);
    //RichtmyerMeshkov test;

    auto eos = test.get_eos({0.0, 0.0, 0.0});

    // Использовать подсеточную реконструкцию начальных данных?
    bool simple_init = false;

    // Распределение консервативных величин
    auto density      = [&test](const Vector3d& r) -> double { return test.density(r); };
    auto momentum_x   = [&test](const Vector3d& r) -> double { return test.density(r) * test.velocity(r).x(); };
    auto momentum_y   = [&test](const Vector3d& r) -> double { return test.density(r) * test.velocity(r).y(); };
    auto momentum_z   = [&test](const Vector3d& r) -> double { return test.density(r) * test.velocity(r).z(); };
    auto total_energy = [&test](const Vector3d& r) -> double {
        Vector3d v = test.velocity(r);
        return test.density(r) * (test.energy(r) + 0.5 * v.dot(v));
    };

    // Задание начальных данных
    auto init_cells = [&](Mesh& mesh) {
        mesh.for_each([&](Cell& cell) {
            if (simple_init ||
               (cell.const_function(density) &&
                cell.const_function(total_energy))) {

                Vector3d r = cell.center();
                cell(U).density  = test.density(r);
                cell(U).velocity = test.velocity(r);
                cell(U).energy   = test.energy(r);
                cell(U).pressure = test.pressure(r);
            }
            else {
                int n = 20;
                double V = cell.volume();

                QState q(
                        cell.integrate_low(density, n) / V, {
                                cell.integrate_low(momentum_x, n) / V,
                                cell.integrate_low(momentum_y, n) / V,
                                cell.integrate_low(momentum_z, n) / V
                        },
                        cell.integrate_low(total_energy, n) / V);

                PState z(q, *eos);

                cell(U).set_state(z);
            }
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

    // Генератор сетки (с граничными условиями) дает тест,
    // число ячеек можно задать
    Rectangle gen = test.generator;
    gen.set_nx(50);

    // Создать сетку
    EuMesh mesh(U, &gen);

    // Создать решатель
    SmFluid solver(eos);
    solver.set_accuracy(2);
    solver.set_CFL(0.5);
    solver.set_limiter("MC");
    solver.set_method(Fluxes::HLLC_M);

    // Сеточная адаптация
    mesh.set_max_level(3);
    mesh.set_distributor(solver.distributor());

    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        init_cells(mesh);
        solver.set_flags(mesh);
        mesh.refine();
    }
    init_cells(mesh);

    size_t n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    while (curr_time < test.max_time()) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";

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

    // Сохранить данные как текст
    CsvFile csv("test2D.csv", 5, pvd.variables);
    csv.save(mesh);

    return 0;
}
