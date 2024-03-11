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

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::mmf;

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
    auto matL = IdealGas::create(1.5, 718.0_J_kgK);
    auto matR = IdealGas::create(2, 718.0_J_kgK);

    Materials mixture;
    mixture += matL;
    mixture += matR;

    // density, velocity, pressure, energy, temperature, mass_frac
    PState zL(0.930887, {-1.98519, 0, 0}, 0.372293, 0.799867, 0.00111402, {1, 0, 0, 0, 0});
    PState zR(0.930951, {1.98507, 0, 0}, 0.372315, 0.399941, 0.00055702, {5.27393e-05, 0.999947, 0, 0, 0});

    auto m_nf = NumFlux::create(Fluxes::GODUNOV);
    auto loc_flux = m_nf->mm_flux(zL, zR, mixture);
    std::cout << loc_flux;
}

int main() {
    threads::on();
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

//    RiemannTesterWithSolver(Fluxes::GODUNOV, test1, 100, 2);
//    RiemannTesterWithSolver(Fluxes::GODUNOV, test3, 100, 2);
//    RiemannTesterWithSolver(Fluxes::GODUNOV, test4, 100, 2);
//    RiemannTesterWithSolver(Fluxes::GODUNOV, test5, 100, 2);
//    RiemannTesterWithSolver(Fluxes::GODUNOV, mm_toro, 100, 2);
//    RiemannTesterWithSolver(Fluxes::GODUNOV, mm_sod, 100, 2);

//    Stopwatch solve;
//    solve.start();
//    RiemannTesterWithSolver2D(Fluxes::GODUNOV, 20, 2, "output_2D_adaptive_2");
//    solve.stop();
//    std::cout << "Time: " << solve.milliseconds();
//    RiemannTesterWithSolver2D(Fluxes::GODUNOV, 20, 1, "output_2D_1");

    threads::off();

    return 0;
}
