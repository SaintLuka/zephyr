/// @file smf_test_axial.cpp
/// @brief Сферический взрыв на двумерной адаптивной сетке с осевой симметрией.

#include <iostream>
#include <iomanip>

#include <zephyr/mesh/mesh.h>
#include <zephyr/geom/generator/rectangle.h>

#include <zephyr/phys/tests/sedov_blast.h>

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
using zephyr::math::SmFluid;
using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;
using zephyr::geom::generator::Rectangle;


// Для быстрого доступа по типу
SmFluid::State U;

/// Переменные для сохранения
double get_lvl(AmrStorage::Item& cell) { return cell.level; }
double get_rho(AmrStorage::Item& cell) { return cell(U).density; }
double get_u(AmrStorage::Item& cell) { return cell(U).velocity.x(); }
double get_v(AmrStorage::Item& cell) { return cell(U).velocity.y(); }
double get_p(AmrStorage::Item& cell) { return cell(U).pressure; }
double get_e(AmrStorage::Item& cell) { return cell(U).energy; }
double get_volume(const AmrStorage::Item &cell) { return cell.get_volume(); }

// Критерий адаптации подобран под задачу.
// Адаптация ячеек с плотностью выше 1.5.
void set_flags(Mesh &mesh) {
    if (!mesh.is_adaptive()) return;

    mesh.for_each([](Cell& cell) {
        const double threshold = 1.5;

        if (cell(U).density > threshold) {
            cell.set_flag(1);
            return;
        }
        for (auto face: cell.faces()) {
            if (face.neib(U).density > threshold) {
                cell.set_flag(1);
                return;
            }
        }
        cell.set_flag(-1);
    });
}

// Массив центров ячеек и граней
NodeStorage centers(EuMesh& mesh) {
    size_t count = 0;
    for (auto& cell: mesh) {
        ++count;
        for (auto& face: cell.faces()) {
            ++count;
        }
    }

    NodeStorage nodes_plain(count, 0);

    count = 0;
    for (auto& cell: mesh) {
        nodes_plain[count].coords = cell.center();
        ++count;
        for (auto& face: cell.faces()) {
            nodes_plain[count].coords = face.center();
            ++count;
        }
    }
    return nodes_plain;
}

int main() {
    mpi::handler init;
    threads::on();

    // Тестовая задача
    SedovBlast3D test({.gamma=1.4, .rho0=1.0, .E=1.0});

    // Уравнение состояния
    IdealGas::Ptr eos = test.eos;

    // Задание начальных данных
    double t0 = test.time_by_radius(0.4);
    test.finish = test.time_by_radius(0.9);

    auto init_cells = [&](Mesh& mesh) {
        mesh.for_each([&](Cell& cell) {
            double r = cell.center().norm();
            Vector3d n = cell.center() / r;
            cell(U).density  = test.density (r, t0);
            cell(U).velocity = test.velocity(r, t0) * n;
            cell(U).pressure = std::max(test.pressure(r, t0), 1.0e-3);
            cell(U).energy   = eos->energy_rP(cell(U).density, cell(U).pressure);
        });
    };

    // Файл для записи
    PvdFile pvd("Sedov", "output");
    pvd.unique_nodes = false;

    size_t n_step = 0;
    double curr_time = t0;
    double next_write = t0;

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"v", get_v};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"volume_as", get_volume};

    pvd.variables += {"rho_exact",
                      [&test, &curr_time](const AmrStorage::Item &cell) -> double {
                          return test.density(cell.center.norm(), curr_time);
                      }};
    pvd.variables += {"v_exact",
                      [&test, &curr_time](const AmrStorage::Item &cell) -> double {
                          return test.velocity(cell.center.norm(), curr_time);
                      }};
    pvd.variables += {"p_exact",
                      [&test, &curr_time](const AmrStorage::Item &cell) -> double {
                          return test.pressure(cell.center.norm(), curr_time);
                      }};
    pvd.variables += {"e_exact",
                      [&test, &curr_time](const AmrStorage::Item &cell) -> double {
                          return test.energy(cell.center.norm(), curr_time);
                      }};

    // Генератор сетки (с граничными условиями) дает тест,
    // число ячеек можно задать
    Rectangle gen(0.0, 1.0, 0.0, 1.0);
    gen.set_boundaries({.left=Boundary::WALL, .right=Boundary::ZOE,
                        .bottom=Boundary::WALL, .top=Boundary::ZOE});
    gen.set_nx(10);
    gen.set_axial(true);

    // Создать сетку
    EuMesh mesh(gen, U);
    mesh.set_decomposition("XY");

    // Создать решатель
    SmFluid solver(eos);
    solver.set_axial(true);
    solver.set_accuracy(2);
    solver.set_CFL(0.5);
    solver.set_limiter("MC");
    solver.set_method(Fluxes::HLLC_M);

    // Сеточная адаптация
    mesh.set_max_level(4);
    mesh.set_distributor(solver.distributor());

    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        init_cells(mesh);
        set_flags(mesh);
        mesh.refine();
    }
    init_cells(mesh);

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
        set_flags(mesh);
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
