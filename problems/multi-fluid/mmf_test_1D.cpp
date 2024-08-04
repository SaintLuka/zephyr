/// @file Одномерные задачи газодинамики

#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/strip.h>
#include <zephyr/mesh/mesh.h>

#include <zephyr/phys/tests/toro.h>
#include <zephyr/phys/tests/rarefied_water.h>

#include <zephyr/math/solver/riemann.h>
#include <zephyr/math/solver/mm_fluid.h>

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
MmFluid::State U;

/// Переменные для сохранения
double get_rho(AmrStorage::Item& cell) { return cell(U).rho; }
double get_u(AmrStorage::Item& cell) { return cell(U).v.x(); }
double get_p(AmrStorage::Item& cell) { return cell(U).p; }
double get_e(AmrStorage::Item& cell) { return cell(U).e; }
double get_T(AmrStorage::Item &cell) { return cell(U).t; }
double get_frac1(AmrStorage::Item &cell) { return cell(U).mass_frac[0]; }
double get_frac2(AmrStorage::Item &cell) { return cell(U).mass_frac[1]; }
double get_c1(AmrStorage::Item &cell) { return cell(U).speeds[0]; }
double get_c2(AmrStorage::Item &cell) { return cell(U).speeds[1]; }
double get_rho1(AmrStorage::Item &cell) { return cell(U).densities[0]; }
double get_rho2(AmrStorage::Item &cell) { return cell(U).densities[1]; }


int main() {
    // Тестовая задача
    //RarefiedWater test;
    ToroTest test(4, true);

    // Состояния слева и справа в тесте
    double Ox = 100.0;
    PState zL(test.density(-Ox), test.velocity(-Ox),
              test.pressure(-Ox), test.energy(-Ox));

    PState zR(test.density(Ox), test.velocity(Ox),
              test.pressure(Ox), test.energy(Ox));

    // Вытащим УрС
    Eos::Ptr eos_L = test.get_eos(-Ox);
    Eos::Ptr eos_R = test.get_eos(+Ox);

    Materials mixture;
    mixture += eos_L;
    mixture += eos_R;

    // Точное решение задачи Римана
    StiffenedGas sgL = eos_L->stiffened_gas(zL.density, zL.pressure);
    StiffenedGas sgR = eos_R->stiffened_gas(zR.density, zR.pressure);
    RiemannSolver exact(zL, zR, sgL, sgR, test.get_x_jump());

    // Файл для записи
    PvdFile pvd("test1D", "output");

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"T", get_T};
    pvd.variables += {"b1", get_frac1};
    pvd.variables += {"b2", get_frac2};
    pvd.variables += {"rho1", get_rho1};
    pvd.variables += {"rho2", get_rho2};

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
    pvd.variables += {"c_exact",
                      [&exact, &curr_time](const AmrStorage::Item &cell) -> double {
                          return exact.sound_speed(cell.center.x(), curr_time);
                      }};
    pvd.variables += {"b_exact",
                      [&exact, &curr_time](const AmrStorage::Item &cell) -> double {
                          return exact.fraction(cell.center.x(), curr_time);
                      }};

    // Создаем одномерную сетку
    Strip gen(test.xmin(), test.xmax());
    gen.set_size(500);

    // Создать сетку
    EuMesh mesh(U, &gen);

    // Создать решатель
    MmFluid solver(mixture, Fluxes::HLL);
    solver.set_CFL(0.9);
    solver.set_acc(1);

    for (auto cell: mesh) {
        cell(U).rho = test.density(cell.center());
        cell(U).v = test.velocity(cell.center());
        cell(U).p = test.pressure(cell.center());
        cell(U).e = test.energy(cell.center());

        cell(U).mass_frac[0] = cell.x() < test.x_jump ? 1.0 : 0.0;
        cell(U).mass_frac[1] = cell.x() < test.x_jump ? 0.0 : 1.0;

        cell(U).t = test.get_eos(cell.x())->temperature_rp(cell(U).rho, cell(U).p);
    }

    double next_write = 0.0;
    size_t n_step = 0;

    while (curr_time <= 1.01 * test.max_time()) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);
            next_write += test.max_time() / 200;
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
    double r_err = 0.0, u_err = 0.0, p_err = 0.0, e_err = 0.0;
    double r_avg = 0.0, u_avg = 0.0, p_avg = 0.0, e_avg = 0.0;
    for (auto cell: mesh) {
        double x = cell.center().x();
        double V = cell.volume();

        r_err += V * std::abs(cell(U).rho   - exact.density(x, curr_time));
        u_err += V * std::abs(cell(U).v.x() - exact.velocity(x, curr_time));
        p_err += V * std::abs(cell(U).p     - exact.pressure(x, curr_time));
        e_err += V * std::abs(cell(U).e     - exact.energy(x, curr_time));

        r_avg += V * std::abs(cell(U).rho  );
        u_avg += V * std::abs(cell(U).v.x());
        p_avg += V * std::abs(cell(U).p    );
        e_avg += V * std::abs(cell(U).e    );
    }
    r_err /= r_avg;
    u_err /= u_avg;
    p_err /= p_avg;
    e_err /= e_avg;

    std::cout << "\nMean errors\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "    Density:     " << r_err << "\n";
    std::cout << "    Velocity:    " << u_err << "\n";
    std::cout << "    Pressure:    " << p_err << "\n";
    std::cout << "    Energy:      " << e_err << "\n";

    return 0;
}
