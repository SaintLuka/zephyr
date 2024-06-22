/// @file Одномерные задачи газодинамики

#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/strip.h>
#include <zephyr/mesh/mesh.h>

#include <zephyr/phys/tests/sod.h>
#include <zephyr/phys/tests/toro.h>
#include <zephyr/phys/tests/quirck.h>
#include <zephyr/phys/tests/shu_osher.h>

#include <zephyr/math/solver/riemann.h>
#include <zephyr/math/solver/sm_fluid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/csv_file.h>

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

using zephyr::mesh::generator::Strip;
using zephyr::math::RiemannSolver;
using zephyr::mesh::EuMesh;


// Для быстрого доступа по типу
SmFluid::State U;

/// Переменные для сохранения
double get_rho(AmrStorage::Item& cell) { return cell(U).rho; }
double get_u(AmrStorage::Item& cell) { return cell(U).v.x(); }
double get_p(AmrStorage::Item& cell) { return cell(U).p; }
double get_e(AmrStorage::Item& cell) { return cell(U).e; }


int main() {
    // Тестовая задача
    //SodTest test;
    //ToroTest test(2);
    //QuirckTest test;
    ShuOsherTest test;

    // Состояния слева и справа в тесте
    double Ox = 100.0;
    PState zL(test.density(-Ox), test.velocity(-Ox),
              test.pressure(-Ox), test.energy(-Ox));

    PState zR(test.density(Ox), test.velocity(Ox),
              test.pressure(Ox), test.energy(Ox));

    // Точное решение задачи Римана
    // (если есть, не для всех одномерных тестов)
    auto& eos = test.get_eos();
    StiffenedGas sg = eos.stiffened_gas(zL.density, zL.pressure);
    RiemannSolver exact(zL, zR, sg, test.get_x_jump());

    // Файл для записи
    PvdFile pvd("test1D", "output");

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};

    double curr_time = 0.0;

    pvd.variables += {"rho_exact",
                      [&exact, &curr_time](const AmrStorage::Item &cell) -> double {
                          return exact.density(cell.center.x(), curr_time);
                      }};
    pvd.variables += {"u_exact",
                      [&exact, &curr_time](const AmrStorage::Item &cell) -> double {
                          return exact.velocity(cell.center.x(), curr_time);
                      }};
    pvd.variables += {"p_exact",
                      [&exact, &curr_time](const AmrStorage::Item &cell) -> double {
                          return exact.pressure(cell.center.x(), curr_time);
                      }};
    pvd.variables += {"e_exact",
                      [&exact, &curr_time](const AmrStorage::Item &cell) -> double {
                          return exact.energy(cell.center.x(), curr_time);
                      }};
    pvd.variables += {"c",
                      [&eos](AmrStorage::Item& cell) -> double {
                          return eos.sound_speed_rp(cell(U).rho, cell(U).p);
                      }};
    pvd.variables += {"c_exact",
                      [&exact, &curr_time](const AmrStorage::Item &cell) -> double {
                          return exact.sound_speed(cell.center.x(), curr_time);
                      }};

    // Создаем одномерную сетку
    Strip gen(test.xmin(), test.xmax());
    gen.set_size(2000);

    // Создать сетку
    EuMesh mesh(U, &gen);

    // Создать решатель
    SmFluid solver(eos, Fluxes::HLLC);
    solver.set_CFL(0.9);
    solver.set_accuracy(2);

    solver.init_cells(mesh, test);

    double next_write = 0.0;
    size_t n_step = 0;

    while (curr_time <= 1.01 * test.max_time()) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);
            next_write += test.max_time() / 100;
        }

        // Обновляем слои
        solver.update(mesh);

        n_step += 1;
        curr_time = solver.get_time();
    }

    // Сохранить данные как текст
    CsvFile csv("test1D.csv", 8, pvd.variables);
    csv.save(mesh);

    // Расчёт ошибок
    double r_err = 0.0, u_err = 0.0, p_err = 0.0, e_err = 0.0, c_err = 0.0;
    double r_avg = 0.0, u_avg = 0.0, p_avg = 0.0, e_avg = 0.0, c_avg = 0.0;
    for (auto cell: mesh) {
        double x = cell.center().x();
        double V = cell.volume();

        r_err += V * std::abs(cell(U).rho   - exact.density(x, curr_time));
        u_err += V * std::abs(cell(U).v.x() - exact.velocity(x, curr_time));
        p_err += V * std::abs(cell(U).p     - exact.pressure(x, curr_time));
        e_err += V * std::abs(cell(U).e     - exact.energy(x, curr_time));
        c_err += V * std::abs(eos.sound_speed_rp(cell(U).rho, cell(U).p) -
                              exact.sound_speed(x, curr_time));

        r_avg += V * std::abs(cell(U).rho  );
        u_avg += V * std::abs(cell(U).v.x());
        p_avg += V * std::abs(cell(U).p    );
        e_avg += V * std::abs(cell(U).e    );
        c_avg += V * std::abs(eos.sound_speed_rp(cell(U).rho, cell(U).p));
    }
    r_err /= r_avg;
    u_err /= u_avg;
    p_err /= p_avg;
    e_err /= e_avg;
    c_err /= c_avg;

    std::cout << "\nMean errors\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "    Density:     " << r_err << "\n";
    std::cout << "    Velocity:    " << u_err << "\n";
    std::cout << "    Pressure:    " << p_err << "\n";
    std::cout << "    Energy:      " << e_err << "\n";
    std::cout << "    Sound Speed: " << c_err << "\n";

    return 0;
}
