#include "fast.h"
#include <zephyr/math/solver/mm_fluid.h>

#include <fstream>
#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/tests/sod.h>
#include <zephyr/phys/tests/toro.h>
#include <zephyr/phys/tests/classic_test.h>

#include <zephyr/math/solver/riemann.h>
#include <zephyr/phys/eos/stiffened_gas.h>

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::mmf;

using zephyr::math::RiemannSolver;


MmFluid::State U;

/// Переменные для сохранения
double get_rho(Storage::Item cell) { return cell(U).rho1; }

double get_u(Storage::Item cell) { return cell(U).v1.x(); }

double get_v(Storage::Item cell) { return cell(U).v1.y(); }

double get_w(Storage::Item cell) { return cell(U).v1.z(); }

double get_p(Storage::Item cell) { return cell(U).p1; }

double get_e(Storage::Item cell) { return cell(U).e1; }


std::vector<double>
RiemannTesterWithSolver(const ClassicTest &test, Fluxes flux, const std::string &filename = "output") {
    // Уравнение состояния
    const Eos &eos = test.get_eos();
    // StiffenedGas eos(1.367, 0.113, 0.273);
    StiffenedGas sg = eos.stiffened_gas(1.0, 1.0);
    Fractions mass_frac({1});
    // Состояния слева и справа в тесте
    Vector3d Ox = 100.0 * Vector3d::UnitX();
    PState zL(test.density(-Ox), test.velocity(-Ox),
              test.pressure(-Ox), test.energy(-Ox),
              eos.temperature_rp(test.density((-Ox)), test.pressure(-Ox)), mass_frac);

    PState zR(test.density(Ox), test.velocity(Ox),
              test.pressure(Ox), test.energy(Ox),
              eos.temperature_rp(test.density((Ox)), test.pressure(Ox)), mass_frac);

    // Точное решение задачи Римана
    RiemannSolver exact(zL.to_smf(), zR.to_smf(), sg, test.get_x_jump());

    // Файл для записи
    PvdFile pvd("mesh", filename);

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};

    double time = 0.0;

    pvd.variables += {"rho_exact",
                      [&exact, &time](const Storage::Item &cell) -> double {
                          return exact.density(cell.center().x(), time);
                      }};
    pvd.variables += {"u_exact",
                      [&exact, &time](const Storage::Item &cell) -> double {
                          return exact.velocity(cell.center().x(), time);
                      }};
    pvd.variables += {"p_exact",
                      [&exact, &time](const Storage::Item &cell) -> double {
                          return exact.pressure(cell.center().x(), time);
                      }};
    pvd.variables += {"e_exact",
                      [&exact, &time](const Storage::Item &cell) -> double {
                          return exact.energy(cell.center().x(), time);
                      }};
    pvd.variables += {"c",
                      [&eos](Storage::Item cell) -> double {
                          return eos.sound_speed_rp(cell(U).rho1, cell(U).p1);
                      }};
    pvd.variables += {"c_exact",
                      [&exact, &time](const Storage::Item &cell) -> double {
                          return exact.sound_speed(cell.center().x(), time);
                      }};

    // Создаем одномерную сетку
    double H = 0.05 * (test.xmax() - test.xmin());
    Rectangle rect(test.xmin(), test.xmax(), -H, +H);
    int n_cells = 500;
    rect.set_sizes(n_cells, 1);
    rect.set_boundary_flags(
            FaceFlag::WALL, FaceFlag::WALL,
            FaceFlag::WALL, FaceFlag::WALL);

    // Создать сетку
    Mesh mesh(U, &rect);

    MmFluid solver(eos, flux);

    for (auto cell: mesh) {
        cell(U).rho1 = test.density(cell.center());
        cell(U).v1 = test.velocity(cell.center());
        cell(U).p1 = test.pressure(cell.center());
        cell(U).e1 = eos.energy_rp(cell(U).rho1, cell(U).p1);
        cell(U).t1 = eos.temperature_rp(cell(U).rho1, cell(U).p1);
        cell(U).mass_frac1 = mass_frac;
    }

    // Число Куранта
    double CFL = 0.9;
    solver.set_CFL(CFL);

    double next_write = 0.0;

    while (time <= 1.01 * test.max_time()) {
        if (time >= next_write) {
            pvd.save(mesh, time);
            next_write += test.max_time() / 100;
        }

        // Определяем dt
        solver.compute_dt(mesh);

        // Расчет по некоторой схеме
        solver.fluxes(mesh);

        // Обновляем слои
        solver.update(mesh);

        time = solver.get_time();
    }

    // расчёт ошибок
    double rho_err = 0.0, u_err = 0.0, p_err = 0.0, e_err = 0.0, c_err = 0.0;
    time = test.max_time();
    for (auto cell: mesh) {
        double x = cell.center().x();
        rho_err += abs(cell(U).rho1 - exact.density(x, time));
        u_err += abs(cell(U).v1.x() - exact.velocity(x, time));
        p_err += abs(cell(U).p1 - exact.pressure(x, time));
        e_err += abs(cell(U).e1 - exact.energy(x, time));
        c_err += abs(eos.sound_speed_rp(cell(U).rho1, cell(U).p1) - exact.sound_speed(x, time));
    }
    rho_err /= n_cells;
    u_err /= n_cells;
    p_err /= n_cells;
    e_err /= n_cells;
    c_err /= n_cells;

    auto fprint = [](const std::string &name, double value) {
        std::cout << name << ": " << value << '\n';
    };

    std::cout << "Test: " << test.get_name() << ", Flux: " << solver.get_flux_name() << "\n";
    std::cout << "Mean average errors:\n";
    fprint("\tdensity error     ", rho_err);
    fprint("\tu error           ", u_err);
    fprint("\tpressure error    ", p_err);
    fprint("\tenergy error      ", e_err);
    fprint("\tsound speed error ", c_err);
    std::cout << '\n';

    return {rho_err, u_err, p_err, e_err, c_err};
}

int main() {
    // Тестовая задача
    SodTest sod_test;
    //ToroTest toro_test(2);

    std::vector<Fluxes> fluxes;
//    fluxes.push_back(Fluxes::HLL);
//    fluxes.push_back(Fluxes::HLLC);
//    fluxes.push_back(Fluxes::HLLC2);
//    fluxes.push_back(Fluxes::CIR2);
//    fluxes.push_back(Fluxes::RUSANOV);
    fluxes.push_back(Fluxes::GODUNOV);

    std::vector<std::vector<double>> sod_errors(fluxes.size(), std::vector<double>(5));
    for (int i = 0; i < fluxes.size(); ++i)
        sod_errors[i] = RiemannTesterWithSolver(sod_test, fluxes[i]);

//    for (int i = 0; i < fluxes.size(); ++i)
//        sod_errors[i] = RiemannTester(sod_test, nfs[i]);

//    std::cout << '\n';
//
//    std::vector<std::vector<double>> toro_errors(fluxes.size(), std::vector<double>(5));
//    for (int i = 0; i < fluxes.size(); ++i)
//        toro_errors[i] = RiemannTesterWithSolver(toro_test, fluxes[i]);

    return 0;
}
