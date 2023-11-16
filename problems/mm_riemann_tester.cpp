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

double get_frac1(Storage::Item cell) { return cell(U).mass_frac1[0]; }

double get_frac2(Storage::Item cell) { return cell(U).mass_frac1[1]; }

struct MmTest {
    std::shared_ptr<Eos> matL, matR;
    double x_jump;
    double max_time;
    double rhoL, rhoR;
    double pL, pR;
    double uL, uR;
    double x_min = 0.0, x_max = 1.0;

    MmTest(const std::shared_ptr<Eos> &matL, const std::shared_ptr<Eos> &matR,
           double x_jump, double max_time,
           double rhoL, double rhoR,
           double pL, double pR,
           double uL, double uR,
           double x_min, double x_max) : matL(matL),
                                         matR(matR),
                                         x_jump(x_jump),
                                         max_time(max_time),
                                         rhoL(rhoL),
                                         rhoR(rhoR),
                                         pL(pL),
                                         pR(pR),
                                         uL(uL),
                                         uR(uR),
                                         x_min(x_min),
                                         x_max(x_max) {}
};

std::vector<double>
RiemannTesterWithSolver(Fluxes flux, const MmTest &test, int n_cells = 10, int acc = 1, const std::string &filename = "output") {
    // Уравнение состояния
    Materials mixture;

    mixture += test.matL;
    mixture += test.matR;

    double x_jump = test.x_jump;
    double max_time = test.max_time;
    double rhoL = test.rhoL, rhoR = test.rhoR;
    double pL = test.pL, pR = test.pR;
    double uL = test.uL, uR = test.uR;
    double eL = test.matL->energy_rp(rhoL, pL), eR = test.matR->energy_rp(rhoR, pR);
    double tL = test.matL->temperature_rp(rhoL, pL), tR = test.matR->temperature_rp(rhoR, pR);
    double x_min = 0.0, x_max = 1.0;

    Fractions mass_fracL({1, 0});
    Fractions mass_fracR({0, 1});
    // Состояния слева и справа в тесте
    PState zL(rhoL, Vector3d(uL, 0, 0), pL, eL, tL, mass_fracL);
    PState zR(rhoR, Vector3d(uR, 0, 0), pR, eL, tR, mass_fracR);

    std::cout << "ZL: " << zL << "\n" << "zR: " << zR << "\n";

    // Точное решение задачи Римана
    RiemannSolver exact(zL.to_smf(), zR.to_smf(),
                        mixture.stiffened_gas(rhoL, pL, mass_fracL),
                        mixture.stiffened_gas(rhoR, pR, mass_fracR),
                        x_jump);

    // Файл для записи
    PvdFile pvd("mesh", filename);

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"frac1", get_frac1};
    pvd.variables += {"frac2", get_frac2};

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
                      [&mixture](Storage::Item cell) -> double {
                          return mixture.sound_speed_rp(cell(U).rho1, cell(U).p1, cell(U).mass_frac1);
                      }};
    pvd.variables += {"c_exact",
                      [&exact, &time](const Storage::Item &cell) -> double {
                          return exact.sound_speed(cell.center().x(), time);
                      }};

    // Создаем одномерную сетку
    double H = 0.05 * (x_max - x_min);
    Rectangle rect(x_min, x_max, -H, +H);
    rect.set_sizes(n_cells, 1);
    rect.set_boundary_flags(
            FaceFlag::WALL, FaceFlag::WALL,
            FaceFlag::WALL, FaceFlag::WALL);

    // Создать сетку
    Mesh mesh(U, &rect);

    MmFluid solver(mixture, flux);
    solver.set_acc(acc);

    for (auto cell: mesh) {
        if (cell.center().x() < x_jump) {
            cell(U).rho1 = rhoL;
            cell(U).v1 = Vector3d(uL, 0, 0);
            cell(U).p1 = pL;
            cell(U).e1 = eL;
            cell(U).t1 = tL;
            cell(U).mass_frac1 = mass_fracL;
        } else {
            cell(U).rho1 = rhoR;
            cell(U).v1 = Vector3d(uR, 0, 0);
            cell(U).p1 = pR;
            cell(U).e1 = eR;
            cell(U).t1 = tR;
            cell(U).mass_frac1 = mass_fracR;
        }
    }

    // Число Куранта
    double CFL = 0.4;
    solver.set_CFL(CFL);

    double next_write = 0.0;
    int n_writes = 100;
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
//            std::cout << "progress: " << round(100 * time / max_time) << "%\n";
            pvd.save(mesh, time);
            next_write += max_time / n_writes;
        }
        solver.update(mesh);

        time = solver.get_time();
    }

    // расчёт ошибок
    double rho_err = 0.0, u_err = 0.0, p_err = 0.0, e_err = 0.0, c_err = 0.0;
    for (auto cell: mesh) {
        double x = cell.center().x();
        rho_err += abs(cell(U).rho1 - exact.density(x, max_time)) / exact.density(x, max_time);
        u_err += abs(cell(U).v1.x() - exact.velocity(x, max_time)) / exact.velocity(x, max_time);
        p_err += abs(cell(U).p1 - exact.pressure(x, max_time)) / exact.pressure(x, max_time);
        e_err += abs(cell(U).e1 - exact.energy(x, max_time)) / exact.energy(x, max_time);
        c_err += abs(mixture.sound_speed_rp(cell(U).rho1, cell(U).p1, cell(U).mass_frac1) -
                     exact.sound_speed(x, max_time)) / exact.sound_speed(x, max_time);
//        std::cout << "idx: " << cell.b_idx() << ", " << cell(U).to_pstate() << "\n";
    }
    rho_err /= n_cells;
    u_err /= n_cells;
    p_err /= n_cells;
    e_err /= n_cells;
    c_err /= n_cells;

    auto fprint = [](const std::string &name, double value) {
        std::cout << name << ": " << value << '\n';
    };

    std::cout << "MultiMaterial test, " << "Flux: " << solver.get_flux_name() << "\n";
    std::cout << "Mean relative errors:\n";
    fprint("\tdensity error     ", rho_err);
    fprint("\tu error           ", u_err);
    fprint("\tpressure error    ", p_err);
    fprint("\tenergy error      ", e_err);
    fprint("\tsound speed error ", c_err);
    std::cout << '\n';

    return {rho_err, u_err, p_err, e_err, c_err};
}

std::vector<double>
RiemannTesterWithSolver2D(Fluxes flux, int n_cells = 10, int acc = 1, const std::string &filename = "output") {
    // Уравнение состояния
    auto matL = IdealGas::create(1.4, 718.0_J_kgK);
    auto matR = IdealGas::create(1.5, 718.0_J_kgK);

    Materials mixture;
    mixture += matL;
    mixture += matR;

    double max_time = 0.01;
    double rho_in = 10, rho_out = 1;
    double p_in = 1e5, p_out = 1e4;
    double u_in = 0, u_out = 0;
    double e_in = matL->energy_rp(rho_in, p_in), e_out = matR->energy_rp(rho_out, p_out);
    double t_in = matL->temperature_rp(rho_in, p_in), t_out = matR->temperature_rp(rho_out, p_out);
    double x_min = 0.0, x_max = 1.0;

    Fractions mass_frac_in({1, 0});
    Fractions mass_frac_out({0, 1});
    // Состояния слева и справа в тесте
    PState zL(rho_in, Vector3d(u_in, 0, 0), p_in, e_in, t_in, mass_frac_in);
    PState zR(rho_out, Vector3d(u_out, 0, 0), p_out, e_in, t_out, mass_frac_out);

    std::cout << "ZL: " << zL << "\n" << "zR: " << zR << "\n";

    // Создаем одномерную сетку
    double H = x_max - x_min;
    double r = 0.05 * H;
    Rectangle rect(x_min, x_max, -H / 2, +H / 2);
    rect.set_sizes(n_cells, n_cells);
    rect.set_boundary_flags(
            FaceFlag::WALL, FaceFlag::WALL,
            FaceFlag::WALL, FaceFlag::WALL);

    // Точное решение задачи Римана
    RiemannSolver exact(zL.to_smf(), zR.to_smf(),
                        mixture.stiffened_gas(rho_in, p_in, mass_frac_in),
                        mixture.stiffened_gas(rho_out, p_out, mass_frac_out),
                        r);

    // Файл для записи
    PvdFile pvd("mesh", filename);

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"frac1", get_frac1};
    pvd.variables += {"frac2", get_frac2};

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
                      [&mixture](Storage::Item cell) -> double {
                          return mixture.sound_speed_rp(cell(U).rho1, cell(U).p1, cell(U).mass_frac1);
                      }};
    pvd.variables += {"c_exact",
                      [&exact, &time](const Storage::Item &cell) -> double {
                          return exact.sound_speed(cell.center().x(), time);
                      }};

    // Создать сетку
    Mesh mesh(U, &rect);

    MmFluid solver(mixture, flux);
    solver.set_acc(acc);

    Vector3d center = {(x_min + x_max) / 2, 0, 0};
    for (auto cell: mesh) {
        if ((cell.center() - center).norm() < r) {
            cell(U).rho1 = rho_in;
            cell(U).v1 = Vector3d(u_in, 0, 0);
            cell(U).p1 = p_in;
            cell(U).e1 = e_in;
            cell(U).t1 = t_in;
            cell(U).mass_frac1 = mass_frac_in;
        } else {
            cell(U).rho1 = rho_out;
            cell(U).v1 = Vector3d(u_out, 0, 0);
            cell(U).p1 = p_out;
            cell(U).e1 = e_out;
            cell(U).t1 = t_out;
            cell(U).mass_frac1 = mass_frac_out;
        }
    }

    // Число Куранта
    double CFL = 0.2;
    solver.set_CFL(CFL);

    double next_write = 0.0;
    int n_writes = 100;
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            std::cout << "progress: " << round(100 * time / max_time) << "%\n";
            pvd.save(mesh, time);
            next_write += max_time / n_writes;
        }
        solver.update(mesh);

        time = solver.get_time();
    }

    // расчёт ошибок
    double rho_err = 0.0, u_err = 0.0, p_err = 0.0, e_err = 0.0, c_err = 0.0;
    for (auto cell: mesh) {
        double x = cell.center().x();
        rho_err += abs(cell(U).rho1 - exact.density(x, max_time)) / exact.density(x, max_time);
        u_err += abs(cell(U).v1.x() - exact.velocity(x, max_time)) / exact.velocity(x, max_time);
        p_err += abs(cell(U).p1 - exact.pressure(x, max_time)) / exact.pressure(x, max_time);
        e_err += abs(cell(U).e1 - exact.energy(x, max_time)) / exact.energy(x, max_time);
        c_err += abs(mixture.sound_speed_rp(cell(U).rho1, cell(U).p1, cell(U).mass_frac1) -
                     exact.sound_speed(x, max_time)) / exact.sound_speed(x, max_time);
//        std::cout << "idx: " << cell.b_idx() << ", " << cell(U).to_pstate() << "\n";
    }
    rho_err /= n_cells;
    u_err /= n_cells;
    p_err /= n_cells;
    e_err /= n_cells;
    c_err /= n_cells;

    auto fprint = [](const std::string &name, double value) {
        std::cout << name << ": " << value << '\n';
    };

    std::cout << "MultiMaterial test, " << "Flux: " << solver.get_flux_name() << "\n";
    std::cout << "Mean relative errors:\n";
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
    ToroTest toro_test(3);

    MmTest transfer(IdealGas::create(1.4, 718.0_J_kgK), IdealGas::create(1.5, 718.0_J_kgK),
                    0.2, 0.1, // x_jump, max_time
                    1.0, 1.0, // rho
                    1e4, 1e4, // p
                    10, 10, // u
                    0, 1.0 // x_min, x_max
    );

    MmTest test1(IdealGas::create(1.4, 718.0_J_kgK), IdealGas::create(1.5, 718.0_J_kgK),
                 0.2, 0.1, // x_jump, max_time
                 10.0, 1.0, // rho
                 2e4, 1e4, // p
                 10, 10, // u
                 0, 1.0 // x_min, x_max
    );

    MmTest test2(IdealGas::create(1.4, 718.0_J_kgK), StiffenedGas::create("Water"),
                 0.5, 0.4,
                 10.0, 800.0,
                 1e5, 1e4,
                 10, -2,
                 0, 1.0
    );

    MmTest test3(IdealGas::create(1.4, 718.0_J_kgK), IdealGas::create(1.5, 718.0_J_kgK),
                 0.2, 0.01, // x_jump, max_time
                 1.0, 1.0, // rho
                 2e4, 1e4, // p
                 0, 0, // u
                 0, 1.0 // x_min, x_max
    );

    MmTest test4(IdealGas::create(1.4, 718.0_J_kgK), IdealGas::create(1.5, 718.0_J_kgK),
                 0.4, 0.004, // x_jump, max_time
                 3.0, 1.0, // rho
                 1e4, 2e4, // p
                 100, 0, // u
                 0, 1.0 // x_min, x_max
    );

    MmTest test5(StiffenedGas::create(1.4, 1000, 100, 718.0_J_kgK), StiffenedGas::create(1.5, 1000, 1000, 718.0_J_kgK),
                 0.4, 0.004, // x_jump, max_time
                 3.0, 1.0, // rho
                 1e4, 2e4, // p
                 100, 0, // u
                 0, 1.0 // x_min, x_max
    );

    double x_min = toro_test.xmin(), x_max = toro_test.xmax();
    MmTest mm_toro(IdealGas::create(1.3, 718.0_J_kgK), IdealGas::create(1.7, 718.0_J_kgK),
                   toro_test.x_jump, toro_test.max_time(), // x_jump, max_time
                   toro_test.density(x_min), toro_test.density(x_max), // rho
                   toro_test.pressure(x_min), toro_test.pressure(x_max), // p
                   toro_test.velocity(x_min).x(), toro_test.velocity(x_max).x(), // u
                   x_min, x_max // x_min, x_max
    );

    x_min = sod_test.xmin(), x_max = sod_test.xmax();
    MmTest mm_sod(IdealGas::create(1.4, 718.0_J_kgK), IdealGas::create(1.5, 718.0_J_kgK),
                  sod_test.x_jump, sod_test.max_time(), // x_jump, max_time
                  sod_test.density(x_min), sod_test.density(x_max), // rho
                  sod_test.pressure(x_min), sod_test.pressure(x_max), // p
                  sod_test.velocity(x_min).x(), sod_test.velocity(x_max).x(), // u
                  x_min, x_max // x_min, x_max
    );

    std::vector<Fluxes> fluxes;
//    fluxes.push_back(Fluxes::HLL);
//    fluxes.push_back(Fluxes::HLLC);
//    fluxes.push_back(Fluxes::HLLC2);
//    fluxes.push_back(Fluxes::RUSANOV2);
    fluxes.push_back(Fluxes::GODUNOV);

//    std::vector<std::vector<double>> sod_errors(fluxes.size(), std::vector<double>(5));
//    for (int i = 1; i < 2; i++)
    RiemannTesterWithSolver2D(Fluxes::GODUNOV, 250, 1);

//    for (int i = 0; i < fluxes.size(); ++i)
//        sod_errors[i] = RiemannTester(sod_test, nfs[i]);

//    std::cout << '\n';
//
//    std::vector<std::vector<double>> toro_errors(fluxes.size(), std::vector<double>(5));
//    for (int i = 0; i < fluxes.size(); ++i)
//        toro_errors[i] = RiemannTesterWithSolver(toro_test, fluxes[i]);

    return 0;
}
