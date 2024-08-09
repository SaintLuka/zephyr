/// @file Двумерные газодинамические задачи, которые ставятся в прямоугольной области

#include <iostream>
#include <iomanip>

#include <zephyr/mesh/mesh.h>

#include <zephyr/phys/tests/test_2D.h>
#include <zephyr/phys/tests/sedov.h>

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
    threads::on();

    // Тестовая задача
    Test2D test(6);
    //SedovBlast test;

    // Начальные данные
    auto init_cells = [&test](Mesh& mesh) {
        for (auto cell: mesh) {
            cell(U).density  = test.density(cell.center());
            cell(U).velocity = test.velocity(cell.center());
            cell(U).pressure = test.pressure(cell.center());
            cell(U).energy   = test.energy(cell.center());
        }
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
    gen.set_nx(300);

    // Создать сетку
    EuMesh mesh(U, &gen);

    // Создать решатель
    SmFluid solver(test.get_eos());
    solver.set_accuracy(2);
    solver.set_CFL(0.5);
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
