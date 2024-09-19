/// @file Двумерные газодинамические задачи, которые ставятся в прямоугольной области

#include <iostream>
#include <iomanip>

#include <zephyr/mesh/mesh.h>

#include <zephyr/phys/tests/test_2D.h>
#include <zephyr/phys/tests/sedov.h>
#include <zephyr/phys/tests/toro.h>

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
double get_rho(AmrStorage::Item& cell) { return cell(U).density; }
double get_u(AmrStorage::Item& cell) { return cell(U).velocity.x(); }
double get_v(AmrStorage::Item& cell) { return cell(U).velocity.y(); }
double get_p(AmrStorage::Item& cell) { return cell(U).pressure; }
double get_e(AmrStorage::Item& cell) { return cell(U).energy; }


int main() {
    threads::on(6);

    // Тестовая задача
    //Test2D test(6);
    //SedovBlast test;

    ToroTest2D test(1, 0.3 * M_PI);

    // Начальные данные
    auto init_cells = [&test](Mesh& mesh) {
        mesh.for_each([&](Cell& cell) {
            using zephyr::geom::Quad;
            Quad quad = ((SqQuad&) cell.geom().vertices).reduce();
            double V = cell.volume();

            int n = 20;

            double R = quad.integrate_low([&](const Vector3d& v) -> double {
                return test.density(v);
            }, n) / V;
            double RU = quad.integrate_low([&](const Vector3d& v) -> double {
                return test.density(v) * test.velocity(v).x();
            }, n) / V;
            double RV = quad.integrate_low([&](const Vector3d& v) -> double {
                return test.density(v) * test.velocity(v).y();
            }, n) / V;
            double RW = quad.integrate_low([&](const Vector3d& v) -> double {
                return test.density(v) * test.velocity(v).z();
            }, n) / V;
            double RE = quad.integrate_low([&](const Vector3d& v) -> double {
                return test.density(v) * test.energy(v);
            }, n) / V;

            cell(U).density  = R;
            cell(U).velocity = Vector3d{RU / R, RV / R, RW / R};
            cell(U).energy   = RE / R;
            cell(U).pressure = test.get_eos().pressure_re(R, RE / R);
        });
    };

    // Файл для записи
    PvdFile pvd("test2D", "output");
    pvd.unique_nodes = true;

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"v", get_v};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};

    // Генератор сетки (с граничными условиями)
    // дает тест, число ячеек можно настроить
    Rectangle gen = test.generator();
    gen.set_nx(800);

    // Создать сетку
    EuMesh mesh(U, &gen);

    // Создать решатель
    SmFluid solver(test.get_eos());
    solver.set_accuracy(2);
    solver.set_CFL(0.5);
    solver.set_limiter("Koren");
    solver.set_method(Fluxes::HLLC);

    // Сеточная адаптация
    mesh.set_max_level(0);
    mesh.set_distributor(solver.distributor());
    
    double curr_time = 0.0;
    double next_write = 0.0;
    size_t n_step = 0;

    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        init_cells(mesh);
        solver.set_flags(mesh);
        mesh.refine();
    }
    init_cells(mesh);

    while (curr_time <= 1.01 * test.max_time()) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";

            pvd.save(mesh, curr_time);
            next_write += test.max_time() / 100;
        }

        // Обновляем слои
        solver.update(mesh);
        solver.set_flags(mesh);
        mesh.refine();

        curr_time += solver.dt();
        n_step += 1;
    }

    // Сохранить данные как текст
    CsvFile csv("test2D.csv", 5, pvd.variables);
    csv.save(mesh);

    return 0;
}
