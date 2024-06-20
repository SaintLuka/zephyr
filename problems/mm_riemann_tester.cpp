#include "fast.h"
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/solver/mm_fluid.h>

#include <fstream>
#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/tests/sod.h>
#include <zephyr/phys/tests/toro.h>
#include <zephyr/phys/tests/classic_test.h>

#include <zephyr/math/solver/riemann.h>
#include <zephyr/phys/eos/stiffened_gas.h>
#include <zephyr/phys/eos/mie_gruneisen.h>
#include <zephyr/utils/matplotlib.h>


using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::mmf;
namespace plt = zephyr::utils::matplotlib;

using zephyr::math::RiemannSolver;


MmFluid::State U;

/// Переменные для сохранения
double get_rho(AmrStorage::Item &cell) { return cell(U).rho; }

double get_u(AmrStorage::Item &cell) { return cell(U).v.x(); }

double get_v(AmrStorage::Item &cell) { return cell(U).v.y(); }

double get_w(AmrStorage::Item &cell) { return cell(U).v.z(); }

double get_p(AmrStorage::Item &cell) { return cell(U).p; }

double get_e(AmrStorage::Item &cell) { return cell(U).e; }

double get_T(AmrStorage::Item &cell) { return cell(U).t; }

double get_frac1(AmrStorage::Item &cell) { return cell(U).mass_frac[0]; }

double get_frac2(AmrStorage::Item &cell) { return cell(U).mass_frac[1]; }

double get_c1(AmrStorage::Item &cell) { return cell(U).speeds[0]; }

double get_c2(AmrStorage::Item &cell) { return cell(U).speeds[1]; }

double get_rho1(AmrStorage::Item &cell) { return cell(U).densities[0]; }

double get_rho2(AmrStorage::Item &cell) { return cell(U).densities[1]; }

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
    double x_min = test.x_min, x_max = test.x_max;

    Fractions mass_fracL({1, 0});
    Fractions mass_fracR({0, 1});
    // Состояния слева и справа в тесте
    PState zL(rhoL, Vector3d(uL, 0, 0), pL, eL, tL, mass_fracL);
    PState zR(rhoR, Vector3d(uR, 0, 0), pR, eR, tR, mass_fracR);

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
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.density(cell.center.x(), time);
                      }};
    pvd.variables += {"u_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.velocity(cell.center.x(), time);
                      }};
    pvd.variables += {"p_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.pressure(cell.center.x(), time);
                      }};
    pvd.variables += {"e_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.energy(cell.center.x(), time);
                      }};
    pvd.variables += {"c",
                      [&mixture](AmrStorage::Item &cell) -> double {
                          return mixture.sound_speed_rp(cell(U).rho, cell(U).p, cell(U).mass_frac);
                      }};
    pvd.variables += {"c_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.sound_speed(cell.center.x(), time);
                      }};

    // Создаем одномерную сетку
    Strip gen(x_min, x_max);
    gen.set_size(n_cells);
    gen.set_boundaries({.left   = Boundary::ZOE, .right = Boundary::ZOE});

    // Создать сетку
    Mesh mesh(U, &gen);

    MmFluid solver(mixture, flux);
    solver.set_acc(acc);
    solver.set_CFL(0.5);
    solver.set_dim(1);

    for (auto cell: mesh) {
        if (cell.center().x() < x_jump) {
            cell(U).rho = rhoL;
            cell(U).v = Vector3d(uL, 0, 0);
            cell(U).p = pL;
            cell(U).e = eL;
            cell(U).t = tL;
            cell(U).mass_frac = mass_fracL;
        } else {
            cell(U).rho = rhoR;
            cell(U).v = Vector3d(uR, 0, 0);
            cell(U).p = pR;
            cell(U).e = eR;
            cell(U).t = tR;
            cell(U).mass_frac = mass_fracR;
        }
    }

    double next_write = 0.0;
    int n_writes = 500;
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

    std::pair<double, double> rho_err = {0, 0}, u_err = {0, 0}, p_err = {0, 0}, e_err = {0, 0}, c_err = {0, 0}; // {mean_err, relative_err}
    /*
    auto sum_err = [](std::pair<double, double> &err, double pred, double real) -> void {
        err.first += abs(pred - real);
        if (real != 0.0)
            err.second += abs(pred - real) / abs(real);
    };
    for (auto cell: mesh) {
        double x = cell.center().x();
        sum_err(rho_err, cell(U).rho, exact.density(x, max_time));
        sum_err(u_err, cell(U).v.x(), exact.velocity(x, max_time));
        sum_err(p_err, cell(U).p, exact.pressure(x, max_time));
        sum_err(e_err, cell(U).e, exact.energy(x, max_time));
        sum_err(c_err, mixture.sound_speed_rp(cell(U).rho, cell(U).p, cell(U).mass_frac), exact.sound_speed(x, max_time));
    }
    rho_err.first /= n_cells;
    rho_err.second /= n_cells;
    u_err.first /= n_cells;
    u_err.second /= n_cells;
    p_err.first /= n_cells;
    p_err.second /= n_cells;
    e_err.first /= n_cells;
    e_err.second /= n_cells;
    c_err.first /= n_cells;
    c_err.second /= n_cells;

    auto fprint = [](const std::string &name, const std::pair<double, double> &value) {
        std::cout << std::fixed << std::setprecision(3) << name << ": " << value.first << " | " << value.second << '\n';
    };

    std::cout << "MultiMaterial test, " << "Flux: " << solver.get_flux_name() << "\n";
    std::cout << "Mean average errors | relative errors:\n";
    fprint("\tdensity error     ", rho_err);
    fprint("\tu error           ", u_err);
    fprint("\tpressure error    ", p_err);
    fprint("\tenergy error      ", e_err);
    fprint("\tsound speed error ", c_err);
    std::cout << '\n';
    */
    return {rho_err.first, u_err.first, p_err.first, e_err.first, c_err.first};
}

void ImpactProblem(int n_cells = 500, int acc = 2, const std::string &filename = "output") {
    // Уравнение состояния
    auto pb = MieGruneisen::create("Pb");
    auto air = IdealGas::create("Air");

    Materials mixture;
    mixture += pb;
    mixture += air;

    double x_min = 0.0_cm, x_max = 1.0_cm;
    double x_jump = 0.5_cm;

    double max_time = 20.0_us;
    double rhoL = 11.4_g_cm3, rhoR = 1.0_kg_m3;
    double pL = 1.0_bar, pR = 1.0_bar;
    double uL = -1200.0_m_s, uR = -1200.0_m_s;
    double eL = pb->energy_rp(rhoL, pL), eR = air->energy_rp(rhoR, pR);
    double tL = pb->temperature_rp(rhoL, pL), tR = air->temperature_rp(rhoR, pR);

    Fractions mass_fracL({1, 0});
    Fractions mass_fracR({0, 1});
    // Состояния слева и справа в тесте
    PState zL(rhoL, Vector3d(uL, 0, 0), pL, eL, tL, mass_fracL);
    PState zR(rhoR, Vector3d(uR, 0, 0), pR, eR, tR, mass_fracR);

    std::cout << "ZL: " << zL << "\n" << "zR: " << zR << "\n";

    // Файл для записи
    PvdFile pvd("mesh", filename);

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"frac1", get_frac1};
    pvd.variables += {"frac2", get_frac2};
    pvd.variables += {"T", get_T};
//    pvd.variables += {"P_min", [&mixture](AmrStorage::Item &cell) -> double {
//        return -mixture.stiffened_gas(cell(U).rho, cell(U).p, cell(U).mass_frac, {.T0 = cell(U).t}).P0;
//    }};
//    pvd.variables += {"gamma", [&mixture](AmrStorage::Item &cell) -> double {
//        return mixture.stiffened_gas(cell(U).rho, cell(U).p, cell(U).mass_frac, {.T0 = cell(U).t}).gamma;
//    }};
//    pvd.variables += {"eps_0", [&mixture](AmrStorage::Item &cell) -> double {
//        return mixture.stiffened_gas(cell(U).rho, cell(U).p, cell(U).mass_frac, {.T0 = cell(U).t}).eps_0;
//    }};

    double time = 0.0;

    // Создаем одномерную сетку
    Strip gen(x_min, x_max);
    gen.set_size(n_cells);
    gen.set_boundaries({.left   = Boundary::WALL, .right = Boundary::ZOE});

    // Создать сетку
    Mesh mesh(U, &gen);

    MmFluid solver(mixture, Fluxes::GODUNOV);
    solver.set_acc(acc);
    solver.set_dim(1);
    solver.set_CFL(0.5);

    for (auto cell: mesh) {
        if (cell.center().x() < x_jump) {
            cell(U).set_state(zL);
        } else {
            cell(U).set_state(zR);
        }
    }

    double next_write = 0.0;
    int n_writes = 2000;
    bool is_set = false;
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            std::cout << "progress: " << std::fixed << std::setprecision(2) << 100 * time / max_time << "%\n";
            pvd.save(mesh, time);
            next_write += max_time / n_writes;
        }

        if (!is_set && time >= 3.6e-6) {
            solver.set_CFL(0.5);
            is_set = true;
        }

        solver.update(mesh);

        time = solver.get_time();
    }
}


std::vector<double>
RiemannTesterWithSolverCSV(Fluxes flux, const MmTest &test, int n_cells = 10, int acc = 1, const std::string &filename = "output.csv") {
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
    PState zR(rhoR, Vector3d(uR, 0, 0), pR, eR, tR, mass_fracR);

    std::cout << "ZL: " << zL << "\n" << "zR: " << zR << "\n";

    // Точное решение задачи Римана
    RiemannSolver exact(zL.to_smf(), zR.to_smf(),
                        mixture.stiffened_gas(rhoL, pL, mass_fracL),
                        mixture.stiffened_gas(rhoR, pR, mass_fracR),
                        x_jump);

    // Файл для записи
    PvdFile pvd("mesh", "1D_tests");

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"frac1", get_frac1};
    pvd.variables += {"frac2", get_frac2};

    double time = 0.0;

    pvd.variables += {"rho_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.density(cell.center.x(), time);
                      }};
    pvd.variables += {"u_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.velocity(cell.center.x(), time);
                      }};
    pvd.variables += {"p_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.pressure(cell.center.x(), time);
                      }};
    pvd.variables += {"e_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.energy(cell.center.x(), time);
                      }};
    pvd.variables += {"c",
                      [&mixture](AmrStorage::Item &cell) -> double {
                          return mixture.sound_speed_rp(cell(U).rho, cell(U).p, cell(U).mass_frac);
                      }};
    pvd.variables += {"c_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.sound_speed(cell.center.x(), time);
                      }};

    // Создаем одномерную сетку
    Strip gen(x_min, x_max);
    gen.set_size(n_cells);
    gen.set_boundaries({.left   = Boundary::ZOE, .right = Boundary::ZOE});

    // Создать сетку
    Mesh mesh(U, &gen);

    MmFluid solver(mixture, flux);
    solver.set_CFL(0.1);
    solver.set_acc(acc);

    for (auto cell: mesh) {
        if (cell.center().x() < x_jump) {
            cell(U).set_state(zL);
        } else {
            cell(U).set_state(zR);
        }
    }

    double next_write = 0.0;
    int n_writes = 100;
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            next_write += max_time / n_writes;
        }
        solver.update(mesh);

        time = solver.get_time();
    }

    CsvFile::save(R"(C:\Users\Diablo\CLionProjects\zephyr\python\output\)" + filename, mesh, 5, pvd.variables);

    // расчёт ошибок
    std::pair<double, double> rho_err = {0, 0}, u_err = {0, 0}, p_err = {0, 0}, e_err = {0, 0}, c_err = {0, 0}; // {mean_err, relative_err}
    auto sum_err = [](std::pair<double, double> &err, double pred, double real) -> void {
        err.first += abs(pred - real);
        if (real != 0.0)
            err.second += abs(pred - real) / abs(real);
    };
    for (auto cell: mesh) {
        double x = cell.center().x();
        sum_err(rho_err, cell(U).rho, exact.density(x, max_time));
        sum_err(u_err, cell(U).v.x(), exact.velocity(x, max_time));
        sum_err(p_err, cell(U).p, exact.pressure(x, max_time));
        sum_err(e_err, cell(U).e, exact.energy(x, max_time));
        sum_err(c_err, mixture.sound_speed_rp(cell(U).rho, cell(U).p, cell(U).mass_frac), exact.sound_speed(x, max_time));
    }
    rho_err.first /= n_cells;
    rho_err.second /= n_cells;
    u_err.first /= n_cells;
    u_err.second /= n_cells;
    p_err.first /= n_cells;
    p_err.second /= n_cells;
    e_err.first /= n_cells;
    e_err.second /= n_cells;
    c_err.first /= n_cells;
    c_err.second /= n_cells;

    auto fprint = [](const std::string &name, const std::pair<double, double> &value) {
        std::cout << std::fixed << std::setprecision(3) << name << ": " << value.first << " | " << value.second << '\n';
    };

    std::cout << "MultiMaterial test, " << "Flux: " << solver.get_flux_name() << "\n";
    std::cout << "Mean average errors | relative errors:\n";
    fprint("\tdensity error     ", rho_err);
    fprint("\tu error           ", u_err);
    fprint("\tpressure error    ", p_err);
    fprint("\tenergy error      ", e_err);
    fprint("\tsound speed error ", c_err);
    std::cout << '\n';

    return {rho_err.first, u_err.first, p_err.first, e_err.first, c_err.first};
}

void ExactSolutionCSV(Fluxes flux, const MmTest &test, int n_cells = 10, const std::string &filename = "output.csv") {
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
    PState zR(rhoR, Vector3d(uR, 0, 0), pR, eR, tR, mass_fracR);

    std::cout << "ZL: " << zL << "\n" << "zR: " << zR << "\n";

    // Точное решение задачи Римана
    RiemannSolver exact(zL.to_smf(), zR.to_smf(),
                        mixture.stiffened_gas(rhoL, pL, mass_fracL),
                        mixture.stiffened_gas(rhoR, pR, mass_fracR),
                        x_jump);

    // Файл для записи
    PvdFile pvd("mesh", "1D_tests");

    double time = max_time;

    pvd.variables += {"rho_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.density(cell.center.x(), time);
                      }};
    pvd.variables += {"u_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.velocity(cell.center.x(), time);
                      }};
    pvd.variables += {"p_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.pressure(cell.center.x(), time);
                      }};
    pvd.variables += {"e_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.energy(cell.center.x(), time);
                      }};
    pvd.variables += {"c",
                      [&mixture](AmrStorage::Item &cell) -> double {
                          return mixture.sound_speed_rp(cell(U).rho, cell(U).p, cell(U).mass_frac);
                      }};
    pvd.variables += {"c_exact",
                      [&exact, &time](AmrStorage::Item &cell) -> double {
                          return exact.sound_speed(cell.center.x(), time);
                      }};

    // Создаем одномерную сетку
    Strip gen(x_min, x_max);
    gen.set_size(n_cells);
    gen.set_boundaries({.left   = Boundary::ZOE, .right = Boundary::ZOE});

    // Создать сетку
    Mesh mesh(U, &gen);

    for (auto cell: mesh) {
        if (cell.center().x() < x_jump) {
            cell(U).set_state(zL);
        } else {
            cell(U).set_state(zR);
        }
    }

    CsvFile::save(R"(C:\Users\Diablo\CLionProjects\zephyr\python\output\)" + filename, mesh, 5, pvd.variables);
}

void calcTests() {
    auto run_exact_test = [](int num, double g1, double g2) {
        ToroTest toro_test(num);
        double x_min = toro_test.xmin(), x_max = toro_test.xmax();
        // gamma, p_inf, eps_0, cV
        MmTest mm_toro(StiffenedGas::create(g1, 0, 0, 718.0_J_kgK), StiffenedGas::create(g2, 0, 0, 718.0_J_kgK),
                       toro_test.x_jump, toro_test.max_time(), // x_jump, max_time
                       toro_test.density(x_min), toro_test.density(x_max), // rho
                       toro_test.pressure(x_min), toro_test.pressure(x_max), // p
                       toro_test.velocity(x_min).x(), toro_test.velocity(x_max).x(), // u
                       x_min, x_max // x_min, x_max
        );
        std::string filename = "exact_toro_test_" + std::to_string(num) + "_" +
                               std::to_string(g1).substr(0, 3) + "_" + std::to_string(g2).substr(0, 3) + ".csv";
        ExactSolutionCSV(Fluxes::GODUNOV, mm_toro, 1000, filename);
    };

    auto run_test = [](int num, double g1, double g2) {
        ToroTest toro_test(num);
        double x_min = toro_test.xmin(), x_max = toro_test.xmax();
        // gamma, p_inf, eps_0, cV
        MmTest mm_toro(StiffenedGas::create(g1, 0, 0, 718.0_J_kgK), StiffenedGas::create(g2, 0, 0, 718.0_J_kgK),
                       toro_test.x_jump, toro_test.max_time(), // x_jump, max_time
                       toro_test.density(x_min), toro_test.density(x_max), // rho
                       toro_test.pressure(x_min), toro_test.pressure(x_max), // p
                       toro_test.velocity(x_min).x(), toro_test.velocity(x_max).x(), // u
                       x_min, x_max // x_min, x_max
        );
        if (num == 2) {
            mm_toro.uL = -1;
            mm_toro.uR = 1;
        }
        std::string filename = "toro_test_" + std::to_string(num) + "_" +
                               std::to_string(g1).substr(0, 3) + "_" + std::to_string(g2).substr(0, 3) + ".csv";
        RiemannTesterWithSolverCSV(Fluxes::GODUNOV, mm_toro, 1000, 2, filename);
    };

//    for (int i = 2; i <= 3; i++) {
//        ToroTest toro_test(i);
//        double x_min = toro_test.xmin(), x_max = toro_test.xmax();
//        for(double g1 = 1.1; g1 <= 1.51; g1 += 0.2){
//            for(double g2 = g1 + 0.1; g2 < 2.41; g2 += 0.2){
//                run_exact_test(i, g1, g2);
//            }
//        }
//    }

//    run_test(2, 1.5, 2.0);
//    run_test(3, 1.1, 1.2);
//    run_test(3, 1.3, 1.4);

    ToroTest toro_test(2);
    double x_min = toro_test.xmin(), x_max = toro_test.xmax();
    // gamma, p_inf, eps_0, cV
    MmTest mm_toro_2(StiffenedGas::create(1.5, 0, 0, 718.0_J_kgK), StiffenedGas::create(2.0, 0, 0, 718.0_J_kgK),
                     toro_test.x_jump, toro_test.max_time(), // x_jump, max_time
                     toro_test.density(x_min), toro_test.density(x_max), // rho
                     toro_test.pressure(x_min), toro_test.pressure(x_max), // p
                     -1, 1, // u
                     x_min, x_max // x_min, x_max
    );
//    RiemannTesterWithSolverCSV(Fluxes::GODUNOV, mm_toro_2, 200, 2, "toro_test2_2_order.csv");
//    RiemannTesterWithSolverCSV(Fluxes::GODUNOV, mm_toro_2, 200, 1, "toro_test2_1_order.csv");

    SodTest sod_test;
    x_min = sod_test.xmin(), x_max = sod_test.xmax();
    MmTest mm_sod(IdealGas::create(1.4, 718.0_J_kgK), IdealGas::create(1.6, 718.0_J_kgK),
                  sod_test.x_jump, sod_test.max_time(), // x_jump, max_time
                  sod_test.density(x_min), sod_test.density(x_max), // rho
                  sod_test.pressure(x_min), sod_test.pressure(x_max), // p
                  sod_test.velocity(x_min).x(), sod_test.velocity(x_max).x(), // u
                  x_min, x_max // x_min, x_max
    );
    RiemannTesterWithSolverCSV(Fluxes::GODUNOV, mm_sod, 200, 2, "sod_test_2_order.csv");
    RiemannTesterWithSolverCSV(Fluxes::GODUNOV, mm_sod, 200, 1, "sod_test_1_order.csv");
}

void twoCellsFlux() {
    auto matL = StiffenedGas::create("Water");
    auto matR = IdealGas::create("Air");

    Materials mixture;
    mixture += matL;
    mixture += matR;

    // density,  velocity, pressure, energy, temperature, mass_frac
//    PState zL(952.086, {1, 0, 0}, -1.026 * 1e7, 812376, 351.429, {1, 0, 0, 0, 0});
//    PState zR(1, {0, 0, 0}, 1e5, 250000, 348.189, {0, 1, 0, 0, 0});
//    PState zL(844.421, {1249.39, -58.5863, 0}, -3.02777e+08, 814071, 295.246, {1, 0, 0, 0, 0});
//    PState zR(933.798, {1323.63, -73.5489, 0}, 48.9301, 814647, 353.548, {1 - 1.04705e-08, 1.04705e-08, 0, 0, 0});
//    PState zL(778.501, {-142.966, -767.994, 0}, 5.34596, 823425, 354.973, {1 - 1.17486e-08, 1.17486e-08, 0, 0, 0});
//    PState zR(812.511, {-140.59, -776.187, 0}, -3.79503e+08, 818268, 278.312, {1, 0, 0, 0, 0});
    PState zL(1244.18, {617.159, -0.254738, 0}, 1.32954 * 1e9, 938374, 547.107, {1, 0, 0, 0, 0});
    PState zR(1244.18, {617.159, 0.254738, 0}, 1.32954 * 1e9, 938374, 547.107, {1, 0, 0, 0, 0});
    zL.energy = matL->energy_rp(zL.density, zL.pressure);
    zL.temperature = matL->temperature_rp(zL.density, zL.pressure);
    zR.energy = matR->energy_rp(zR.density, zR.pressure);
    zR.temperature = matR->temperature_rp(zR.density, zR.pressure);

    std::cout << "ZL: " << zL << "\n" << "zR: " << zR << "\n";

    auto m_nf = NumFlux::create(Fluxes::GODUNOV);
    auto loc_flux = m_nf->mm_flux(zL, zR, mixture);
    std::cout << loc_flux;
}

void find_PT() {
//    QState q(952.086, {-48.0352, 0, 0}, 7.76441e+08, FractionsFlux(std::vector<double>{952.086, 3.02876e-05, 0, 0, 0}));
/*
Failed to calc PState from QState in fluxes_stage2
QState: mass: 861.669, momentum: {2.49503e+06, -41249.7, 0}, energy: 4.07937e+09, mass_frac: [856.787, 4.88194, 0, 0, 0]
Failed to calc PState from QState in fluxes_stage2
QState: mass: 879.816, momentum: {2.55266e+06, -40768.9, 0}, energy: 4.14476e+09, mass_frac: [874.83, 4.98573, 0, 0, 0]
PState: density: 861.669, velocity: {2895.58, -47.8719, 0}, pressure: nan, temperature: nan, energy: 540935, mass_frac: [0.994334, 0.00566568, 0, 0, 0]
Previous PState:
density: 855.307, velocity: {2863.68, -50.7467, 0}, pressure: 2.07384e+08, temperature: 420.184, energy: 978081, mass_frac: [0.994294, 0.00570604, 0, 0, 0]

Failed to calc PState from QState in fluxes_stage2
QState: mass: 43.3173, momentum: {-75385.5, 154394, 0}, energy: 1.98132e+08, mass_frac: [39.8422, 3.4751, 0, 0, 0]
PState: density: 43.3173, velocity: {-1740.31, 3564.25, 0}, pressure: nan, temperature: nan, energy: -3.29232e+06, mass_
frac: [0.919776, 0.0802243, 0, 0, 0]
Previous PState: density: 41.7965, velocity: {-1730.7, 3563, 0}, pressure: 5.39294e+06, temperature: 2191.69, energy: 1.
11809e+07, mass_frac: [0.916691, 0.0833089, 0, 0, 0]

Failed to calc PState from QState in fluxes_stage2
QState: mass: 890.801, momentum: {1.19051e+06, -143374, 0}, energy: -1.65055e+09, mass_frac: [890.068, 0.73281, 0, 0, 0]
PState: density: 890.801, velocity: {1336.45, -160.949, 0}, pressure: nan, temperature: nan, energy: -2.75889e+06, mass_
frac: [0.999177, 0.000822642, 0, 0, 0]
Previous PState: density: 885.444, velocity: {1332.07, -161.556, 0}, pressure: 179612, temperature: 307.596, energy: 531
182, mass_frac: [0.999094, 0.000905938, 0, 0, 0]

 Failed to calc PState from QState in fluxes
QState: mass: 897.95, momentum: {891630, -112319, 0}, energy: -1.34419e+09, mass_frac: {897.557, 0.393557, 0, 0, 0}
PState: density: 897.95, velocity: {992.962, -125.084, 0}, pressure: nan, temperature: nan, energy: -1.99776e+06, mass_f
rac: {0.999562, 0.000438284, 0, 0, 0}
Previous PState: density: 892.942, velocity: {988.304, -123.779, 0}, pressure: 1.26422e+06, temperature: 358.184, energy
: 841574, mass_frac: {0.999549, 0.000451239, 0, 0, 0}

 */
//    QState q(879.816, {2.55266e+06, -40768.9, 0}, 4.14476e+09, FractionsFlux(std::vector<double>{874.83, 4.98573, 0, 0, 0}));
//    QState q(43.3173, {-75385.5, 154394, 0}, 1.98132e+08, FractionsFlux(std::vector<double>{39.8422, 3.4751, 0, 0, 0}));
//    QState q(890.801, {1.19051e+06, -143374, 0}, -1.65055e+09, FractionsFlux(std::vector<double>{890.068, 0.73281, 0, 0, 0}));
    QState q(904.759, {992.962, -125.084, 0}, -1.34419e+09, FractionsFlux(std::vector<double>{1, 0, 0, 0, 0}));
    Materials mixture_eos;
    mixture_eos += StiffenedGas::create("Water");
    mixture_eos += StiffenedGas::create("Air");
    std::cout << PState(q, mixture_eos, -2, 300) << "\n";

    PState p_state;
    auto mass_frac = Fractions(q.mass_frac);
    p_state.mass_frac = mass_frac;
    p_state.density = q.mass;
//    p_state.velocity = q.momentum / p_state.density;
//    p_state.energy = q.energy / p_state.density - 0.5 * p_state.velocity.squaredNorm();
    p_state.energy = 824908;

    std::vector<StiffenedGas> mixture = {StiffenedGas("Water"), StiffenedGas("Air")};
    int n = mixture.size();

    double eps0;
    std::vector<double> B(n, 0), g(n), P0(n), T0(n);
    for (int i = 0; i < n; i++) {
        if (mass_frac.has(i)) {
            eps0 += mass_frac[i] * mixture[i].eps_0;
            B[i] = mass_frac[i] * mixture[i].Cv;
        }
        g[i] = mixture[i].gamma - 1;
        P0[i] = mixture[i].P0;
        T0[i] = mixture[i].T0;
    }
    double P_min = mixture_eos.min_pressure(mass_frac);
    std::cout << "P_min: " << P_min << '\n';

    P_min = -3e8;
    double P_max = -5e7;
    int size = 10000;
    double step = (P_max - P_min) / size;

    std::vector<double> P(size);
    std::vector<double> f(size);
    int sol_idx = 0;
    for (int idx = 0; idx < size; idx++) {
        double p = idx == 0 ? P_min : P[idx - 1] + step;
        P[idx] = p;
        double sum1 = 0, sum2 = 0, sum3 = 0;
        for (int i = 0; i < n; i++) {
            double tmp = B[i] * (1 + g[i] * P0[i] / (p + P0[i]));
            sum1 += tmp;

            double sub_sum = 0;
            for (int k = 0; k < n; k++) {
                if (k == i)
                    continue;
                sub_sum += B[k] * g[k] * (T0[k] - T0[i]) / (p + P0[k]);
            }
            sum2 += tmp * sub_sum;

            sum3 += B[i] * g[i] / (p + P0[i]);
        }
        sum1 /= p_state.density;
        sum3 *= (p_state.energy - eps0);
        f[idx] = sum1 + sum2 - sum3;

        if (abs(f[idx]) < abs(f[sol_idx])) {
            sol_idx = idx;
        }
    }

    plt::figure_size(14, 8);
    std::map<std::string, std::string> m({{"color", "orange"}});
    std::cout << "P: " << P[sol_idx] << ", f: " << f[sol_idx] << '\n';
    plt::plot(P, f);
//    std::vector<double> sub_p(P.begin() + sol_idx - size / 1000, P.begin() + sol_idx + size / 1000);
//    std::vector<double> sub_f(f.begin() + sol_idx - size / 1000, f.begin() + sol_idx + size / 1000);
//    plt::plot(sub_p, sub_f);
    plt::tight_layout();
    plt::show();

//    p.pressure = mixture.pressure_re(density, energy, mass_frac, {.P0=P0, .T0=T0});
//    p.temperature = mixture.temperature_rp(density, pressure, mass_frac, {.T0=T0});
}

void find_PT2() {
    /*
    Failed to calc PState from QState in fluxes_stage2
    QState: mass: 2304.05, momentum: {-411216, 0, 0}, energy: -2.08431e+08, mass_frac: {2304.05, 8.50948e-05, 0, 0, 0}
    PState: density: 2304.05, velocity: {-178.476, 0, 0}, pressure: nan, temperature: nan, energy: -106390, mass_frac: {1, 3.69328e-08, 0, 0, 0}
    Previous PState: density: 1736.32, velocity: {-64.1133, 0, 0}, pressure: 74.4769, temperature: 2676.13, energy: 365287, mass_frac: {1, 4.50739e-08, 0, 0, 0}

     Failed to calc PState from QState in fluxes_stage2
    QState: mass: 1344.46, momentum: {-1.73939e+06, 0, 0}, energy: -1.19307e+09, mass_frac: {1344.46, 5.43148e-05, 0, 0, 0}
    PState: density: 1344.46, velocity: {-1293.75, 0, 0}, pressure: nan, temperature: nan, energy: -1.72429e+06, mass_frac:{1, 4.03989e-08, 0, 0, 0}
    Previous PState: density: 1117.4, velocity: {698.168, 0, 0}, pressure: 61.6383, temperature: 3077.63, energy: 453402, mass_frac: {1, 5.43048e-08, 0, 0, 0}

     Failed to calc PState from QState in fluxes_stage2
    QState: mass: 1128.72, momentum: {-2.13144e+06, 0, 0}, energy: -1.12074e+09, mass_frac: {1128.72, 6.31637e-05, 0, 0, 0}
    PState: density: 1128.72, velocity: {-1888.36, 0, 0}, pressure: nan, temperature: nan, energy: -2.77588e+06, mass_frac: {1, 5.59605e-08, 0, 0, 0}
    Previous PState: density: 915.281, velocity: {703.643, 0, 0}, pressure: 110.601, temperature: 4780.4, energy: 920560, mass_frac: {1, 7.63444e-08, 0, 0, 0}

     Failed to calc PState from QState in fluxes_stage2
    QState: mass: 2062.58, momentum: {-1.35768e+06, -420176, 0}, energy: -7.80816e+08, mass_frac: {2062.58, 4.42559e-08, 0, 0, 0}
    PState: density: 2062.58, velocity: {-658.241, -203.713, 0}, pressure: nan, temperature: nan, energy: -615952, mass_frac: {1, 2.14565e-11, 0, 0, 0}
    Previous PState: density: 2088.02, velocity: {-659.333, -204.721, 0}, pressure: 128.072, temperature: 17533.1, energy: -411741, mass_frac: {1, 2.14709e-11, 0, 0, 0}

     Failed to calc PState from QState in fluxes_stage1
    QState: mass: 1680.78, momentum: {-931039, 0, 0}, energy: -4.8223e+08, mass_frac: {1680.78, 0.00708158, 0, 0, 0}
    PState: density: 1680.78, velocity: {-553.931, 0, 0}, pressure: nan, temperature: nan, energy: -440328, mass_frac: {0.999996, 4.21326e-06, 0, 0, 0}
    Previous PState: density: 1399.53, velocity: {-457.989, 0, 0}, pressure: 734232, temperature: 24269.7, energy: -59244.6, mass_frac: {0.999995, 5.19436e-06, 0, 0, 0}
     */
//    QState q(2304.05, {-411216, 0, 0}, -2.08431e+08, FractionsFlux(std::vector<double>{2304.05, 8.50948e-05, 0, 0, 0}));
//    double P0 = 74.4769, T0 = 2676.13;
//    QState q(1344.46, {-1.73939e+06, 0, 0}, -1.19307e+09, FractionsFlux(std::vector<double>{1344.46, 5.43148e-05, 0, 0, 0}));
//    double P0 = 61.6383, T0 = 3077.63;
//    QState q(1128.72, {-2.13144e+06, 0, 0}, -1.12074e+09, FractionsFlux(std::vector<double>{1128.72, 6.31637e-05, 0, 0, 0}));
//    double P0 = 110.601, T0 = 4780.4;
//    QState q(2062.58, {-1.35768e+06, -420176, 0}, -7.80816e+08, FractionsFlux(std::vector<double>{2062.58, 4.42559e-08, 0, 0, 0}));
//    double P0 = 128.072, T0 = 17533.1;
    QState q(1680.78, {-931039, 0, 0}, -4.8223e+08, FractionsFlux(std::vector<double>{1680.78, 0.00708158, 0, 0, 0}));
    double P0 = -1.61141e+10, T0 = 24269.7;
    Materials mixture_eos;
    mixture_eos += MieGruneisen::create("Pb");
    mixture_eos += IdealGas::create("Air");

    std::cout << PState(q, mixture_eos, P0, T0) << "\n";
    std::cout << PState(q, mixture_eos, -P0, T0) << "\n";

    std::cout << PState(q, mixture_eos, P0, 50) << "\n";
    std::cout << PState(q, mixture_eos, -P0, 50) << "\n";

    std::cout << PState(q, mixture_eos, P0, 2000) << "\n";
    std::cout << PState(q, mixture_eos, -P0, 2000) << "\n";

    std::cout << PState(q, mixture_eos, 1e5, T0) << "\n";
    std::cout << PState(q, mixture_eos, -1e5, T0) << "\n";

    std::cout << PState(q, mixture_eos, 1e5, 2000) << "\n";
    std::cout << PState(q, mixture_eos, -1e5, 2000) << "\n";

    std::cout << PState(q, mixture_eos, 1e2, 20) << "\n";
    std::cout << PState(q, mixture_eos, -1e2, 20) << "\n";

    std::cout << PState(q, mixture_eos, 1e5, 20) << "\n";
    std::cout << PState(q, mixture_eos, -1e5, 20) << "\n";

    std::cout << PState(q, mixture_eos, 1e2, 2000) << "\n";
    std::cout << PState(q, mixture_eos, -1e2, 2000) << "\n";
}

int main() {
    threads::on(16);
    // Тестовая задача
    SodTest sod_test;
    ToroTest toro_test(3);
    

    MmTest waterAir(StiffenedGas::create("Water"), IdealGas::create("Air"),
                    0, 1e-3, // x_jump, max_time
                    952.086, 1.0, // rho
                    -1e5, 1e5, // p
                    10, 0, // u
                    -1, 1  // x_min, x_max
    );

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
                 0.5, 0.1,
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

//    twoCellsFlux();

//    RiemannTesterWithSolver(Fluxes::GODUNOV, waterAir, 1000, 1, "water_air");
//    RiemannTesterWithSolver(Fluxes::GODUNOV, test3, 100, 2);
//    RiemannTesterWithSolver(Fluxes::GODUNOV, test4, 100, 2);
//    RiemannTesterWithSolver(Fluxes::GODUNOV, test5, 100, 2);
//    RiemannTesterWithSolver(Fluxes::GODUNOV, mm_toro, 100, 2);
//    RiemannTesterWithSolver(Fluxes::GODUNOV, mm_sod, 100, 2);

//    find_PT2();
    ImpactProblem(1000, 2, "impact_problem2");

//    Stopwatch solve;
//    solve.start();
//    RiemannTesterWithSolver2D(Fluxes::GODUNOV, 20, 2, "output_2D_adaptive_2");
//    solve.stop();
//    std::cout << "Time: " << solve.milliseconds();
//    RiemannTesterWithSolver2D(Fluxes::GODUNOV, 20, 1, "output_2D_1");

    threads::off();

    return 0;
}
