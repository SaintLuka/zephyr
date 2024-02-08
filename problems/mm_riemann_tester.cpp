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

void RiemannTesterWithSolver2D(Fluxes flux, int n_cells = 10, int acc = 1, const std::string &filename = "output") {
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
    PState zR(rho_out, Vector3d(u_out, 0, 0), p_out, e_out, t_out, mass_frac_out);

    std::cout << "ZL: " << zL << "\n" << "zR: " << zR << "\n";

    // Создаем одномерную сетку
    double H = x_max - x_min;
    double r = 0.05 * H;
    Rectangle rect(x_min, x_max, -H / 2, +H / 2);
    rect.set_sizes(n_cells, n_cells);
    rect.set_boundaries({
                                .left   = Boundary::ZOE, .right = Boundary::ZOE,
                                .bottom = Boundary::ZOE, .top   = Boundary::ZOE});

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

    // Создать сетку
    Mesh mesh(U, &rect);

    MmFluid solver(mixture, flux);
    solver.set_acc(acc);

    // Число Куранта
    double CFL = 0.4;
    solver.set_CFL(CFL);

    // Настраиваем адаптацию
    mesh.set_max_level(5);
    mesh.set_distributor(solver.distributor());

    Vector3d center = {(x_min + x_max) / 2, 0, 0};
    // Адаптация под начальные данные
    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        for (auto cell: mesh) {
            if ((cell.center() - center).norm() < r) {
                cell(U).rho = rho_in;
                cell(U).v = Vector3d(u_in, 0, 0);
                cell(U).p = p_in;
                cell(U).e = e_in;
                cell(U).t = t_in;
                cell(U).mass_frac = mass_frac_in;
            } else {
                cell(U).rho = rho_out;
                cell(U).v = Vector3d(u_out, 0, 0);
                cell(U).p = p_out;
                cell(U).e = e_out;
                cell(U).t = t_out;
                cell(U).mass_frac = mass_frac_out;
            }
        }
        solver.set_flags(mesh);
        mesh.refine();
    }

    double next_write = 0.0;
    int n_writes = 100;
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            std::cout << "progress: " << int(round(100 * time / max_time)) << "%\n";
            pvd.save(mesh, time);
            next_write += max_time / n_writes;
        }

        // шаг решения
        solver.update(mesh);

        // Установить флаги адаптации
        solver.set_flags(mesh);

        // Адаптировать сетку
        mesh.refine();

        time = solver.get_time();
    }
}

void RiemannTesterWithSolverVertical(double g = 9.81, int acc = 2, const std::string &filename = "output") {
    // Уравнение состояния
    auto mat_up = IdealGas::create(1.4, 718.0_J_kgK);
    auto mat_down = IdealGas::create(1.6, 718.0_J_kgK);

    Materials mixture;
    mixture += mat_up;
    mixture += mat_down;

    double max_time = 2;
    double rho_up = 100, rho_down = 20;
    double p_up = 1e4, p_down = 1e4;
    double v_up = 0, v_down = 0;
    double e_up = mat_up->energy_rp(rho_up, p_up), e_down = mat_down->energy_rp(rho_down, p_down);
    double t_up = mat_up->temperature_rp(rho_up, p_up), t_down = mat_down->temperature_rp(rho_down, p_down);
    double x_min = 0, x_max = 0.1;
    double y_min = 0.0, y_max = 1.0;
    double y_jump = 0.5 * (y_max - y_min);

    Fractions mass_frac_up({1, 0});
    Fractions mass_frac_down({0, 1});
    // Состояния слева и справа в тесте
    PState z_up(rho_up, Vector3d(0, v_up, 0), p_up, e_up, t_up, mass_frac_up);
    PState z_down(rho_down, Vector3d(0, v_down, 0), p_down, e_down, t_down, mass_frac_down);

    std::cout << "Z_u: " << z_up << "\n" << "Z_d: " << z_down << "\n";

    // Создаем одномерную сетку
    double H = y_max - y_min;
    Rectangle rect(x_min, x_max, 0, H);
    rect.set_sizes(4, 20);
    rect.set_boundaries({
                                .left   = Boundary::WALL, .right = Boundary::WALL,
                                .bottom = Boundary::WALL, .top   = Boundary::WALL});

    // Файл для записи
    PvdFile pvd("mesh", filename);

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"v", get_v};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"frac1", get_frac1};
    pvd.variables += {"frac2", get_frac2};

    double time = 0.0;

    // Создать сетку
    Mesh mesh(U, &rect);

    MmFluid solver(mixture, Fluxes::GODUNOV, g);
    solver.set_acc(acc);

    // Число Куранта
    double CFL = 0.4;
    solver.set_CFL(CFL);

    // Настраиваем адаптацию
    mesh.set_max_level(5);
    mesh.set_distributor(solver.distributor());

    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        for (auto cell: mesh) {
            if (cell.center().y() > y_jump) {
                cell(U).set_state(z_up);
            } else {
                cell(U).set_state(z_down);
            }
        }
        solver.set_flags(mesh);
        mesh.refine();
    }

    double L = x_max - x_min;
    for (auto cell: mesh) {
        double x = cell.center().x(), y = cell.center().y();
        if (cell(U).mass_frac[1] == 0 && y - y_jump < 0.05 * H && abs(x - L / 2) < 0.1 * L) {
            bool exist = false;
            for (auto &face: cell.faces()) {
                if (face.is_boundary())
                    continue;
                exist = face.neib()(U).mass_frac[1] > 0;
                if (exist)
                    break;
            }
            if (exist) {
                cell(U).mass_frac[0] = 0.9;
                cell(U).mass_frac[1] = 0.1;
            }
        }
    }

    double next_write = 0.0;
    int n_writes = 200;
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            std::cout << "progress: " << std::fixed << std::setprecision(1) << 100 * time / max_time << "%\n";
            pvd.save(mesh, time);
            next_write += max_time / n_writes;
        }

        // шаг решения
        solver.update(mesh);

        // Установить флаги адаптации
        solver.set_flags(mesh);

        // Адаптировать сетку
        mesh.refine();

        time = solver.get_time();
    }
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

void KelvinHelmholtzInstability(int acc = 2, const std::string &filename = "output") {
    // Уравнение состояния
    auto mat_up = IdealGas::create(1.4, 718.0_J_kgK);
    auto mat_down = IdealGas::create(1.4, 718.0_J_kgK);

    Materials mixture;
    mixture += mat_up;
    mixture += mat_down;

    double max_time = 0.5;
    double rho_up = 11, rho_down = 10;
    double p_up = 1e4, p_down = 1e4;
    double u_up = 3, u_down = -3;
    double e_up = mat_up->energy_rp(rho_up, p_up), e_down = mat_down->energy_rp(rho_down, p_down);
    double t_up = mat_up->temperature_rp(rho_up, p_up), t_down = mat_down->temperature_rp(rho_down, p_down);
    double x_min = 0, x_max = 0.4;
    double y_min = 0.0, y_max = 0.15;
    double L = x_max - x_min, H = y_max - y_min;
    double y_jump = y_min + 0.5 * H;

    Fractions mass_frac_up({1, 0});
    Fractions mass_frac_down({0, 1});
    // Состояния слева и справа в тесте
    PState z_up(rho_up, Vector3d(u_up, 0, 0), p_up, e_up, t_up, mass_frac_up);
    PState z_down(rho_down, Vector3d(u_down, 0, 0), p_down, e_down, t_down, mass_frac_down);

    std::cout << "Z_u: " << z_up << "\n" << "Z_d: " << z_down << "\n";

    Rectangle rect(x_min, x_max, 0, H);
    rect.set_sizes(6, 10);
    rect.set_boundaries({
                                .left   = Boundary::ZOE, .right = Boundary::ZOE,
                                .bottom = Boundary::ZOE, .top   = Boundary::ZOE});

    // Файл для записи
    PvdFile pvd("mesh", filename);

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"v", get_v};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"frac1", get_frac1};
    pvd.variables += {"frac2", get_frac2};

    double time = 0.0;

    // Создать сетку
    Mesh mesh(U, &rect);

    MmFluid solver(mixture, Fluxes::GODUNOV, 1);
    solver.set_acc(acc);

    // Число Куранта
    double CFL = 0.4;
    solver.set_CFL(CFL);

    // Настраиваем адаптацию
    mesh.set_max_level(6);
    mesh.set_distributor(solver.distributor());
    // y_jump + sin(10 * pi * x / L) * H * 0.1

    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        for (auto cell: mesh) {
            if (cell.center().y() > y_jump + 0.005 * H * sin(M_PI * cell.center().x() / L)) {
                cell(U).set_state(z_up);
            } else {
                cell(U).set_state(z_down);
            }
        }
        solver.set_flags(mesh);
        mesh.refine();
    }

//    for (auto cell: mesh) {
//        double x = cell.center().x(), y = cell.center().y();
//        if (cell(U).mass_frac[1] == 0 && y - y_jump < 0.2 * H && abs(x - L / 2) < 0.1 * L) {
//            bool exist = false;
//            for (auto &face: cell.faces()) {
//                if (face.is_boundary())
//                    continue;
//                exist = face.neib()(U).mass_frac[1] > 0;
//                if (exist)
//                    break;
//            }
//            if (exist) {
////                cell(U).mass_frac[0] = 0.2;
////                cell(U).mass_frac[1] = 0.8;
//                cell(U).v.y() = -2;
//            }
//        }
//    }


    double next_write = 0.0;
    int n_writes = 100;
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            std::cout << "progress: " << std::fixed << std::setprecision(1) << 100 * time / max_time << "%\n";
            pvd.save(mesh, time);
            next_write += max_time / n_writes;
        }

        // шаг решения
        solver.update(mesh);

        // Установить флаги адаптации
        solver.set_flags(mesh);

        // Адаптировать сетку
        mesh.refine();

        time = solver.get_time();
    }
}

void Bubble2D(int n_cells = 20, int acc = 2, const std::string &filename = "output") {
    // Уравнение состояния
    auto water = StiffenedGas::create("Water");
    auto air = IdealGas::create("Air");

    Materials mixture;
    mixture += water;
    mixture += air;

    double max_time = 4.5e-6;
    double rho_wave = 1323.65_kg_m3, rho_water = 1000.0_kg_m3, rho_air = 1.0_kg_m3;
    double p_wave = 1.9e4_bar, p_water = 1.0_bar, p_air = 1.0_bar;
    double u_wave = 681.58_m_s, u_water = 0, u_air = 0;
    double e_wave = water->energy_rp(rho_wave, p_wave), e_water = water->energy_rp(rho_water, p_water), e_air = air->energy_rp(rho_air, p_air);
    double t_wave = water->temperature_rp(rho_wave, p_wave), t_water = water->temperature_rp(rho_water, p_water), t_air = air->temperature_rp(rho_air, p_air);

    double x_min = 0.0_cm, x_max = 0.8_cm;
    double y_min = 0.0_cm, y_max = 0.8_cm;
    double r = 0.3_cm, x_bubble = (x_max + x_min) / 2, y_bubble = (y_max + y_min) / 2;
    double x_wave = 0.06_cm;

    Fractions mass_frac_water({1, 0});
    Fractions mass_frac_air({0, 1});

    PState state_air(rho_air, Vector3d{u_air, 0, 0}, p_air, e_air, t_air, mass_frac_air);
    PState state_wave(rho_wave, Vector3d{u_wave, 0, 0}, p_wave, e_wave, t_wave, mass_frac_water);
    PState state_water(rho_water, Vector3d{u_water, 0, 0}, p_water, e_water, t_water, mass_frac_water);

    // Создаем одномерную сетку
    Rectangle rect(x_min, x_max, y_min, y_max);
    rect.set_sizes(n_cells, n_cells);
    rect.set_boundaries({
                                .left   = Boundary::ZOE, .right = Boundary::ZOE,
                                .bottom = Boundary::WALL, .top   = Boundary::WALL});

    // Файл для записи
    PvdFile pvd("mesh", filename);

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"v", get_v};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"frac1", get_frac1};
    pvd.variables += {"frac2", get_frac2};
    pvd.variables += {"c1", get_c1};
    pvd.variables += {"c2", get_c2};
    pvd.variables += {"rho1", get_rho1};
    pvd.variables += {"rho2", get_rho2};

    double time = 0.0;

    // Создать сетку
    Mesh mesh(U, &rect);

    MmFluid solver(mixture, Fluxes::GODUNOV);
    solver.set_acc(acc);

    // Число Куранта
    double CFL = 0.02;
    solver.set_CFL(CFL);

    // Настраиваем адаптацию
    mesh.set_max_level(5);
    mesh.set_distributor(solver.distributor());

    Vector3d bubble_center = {x_bubble, y_bubble, 0};
    // Адаптация под начальные данные
    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        for (auto cell: mesh) {
            if (cell.center().x() < x_wave) {
                cell(U).set_state(state_wave);
                continue;
            }

            if ((cell.center() - bubble_center).norm() < r) {
                cell(U).set_state(state_air);
            } else {
                cell(U).set_state(state_water);
            }
        }
        solver.set_flags(mesh);
        mesh.refine();
    }

    double next_write = 0.0;
    int n_writes = 200;
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            std::cout << "progress: " << std::fixed << std::setprecision(1) << 100 * time / max_time << "%\n";
            pvd.save(mesh, time);
            next_write += max_time / n_writes;
        }

        // шаг решения
        solver.update(mesh);

        // Установить флаги адаптации
        solver.set_flags(mesh);

        // Адаптировать сетку
        mesh.refine();

        time = solver.get_time();
    }
}

void Bubble2DStatic(int n_cells = 100, int acc = 2, const std::string &filename = "output") {
    // Уравнение состояния
    auto water = StiffenedGas::create("Water");
    auto air = IdealGas::create("Air");

    Materials mixture;
    mixture += water;
    mixture += air;
//    Failed to calc PState from QState
//    QState: mass: 910.166, momentum: {59551.3, -312854, 0}, energy: 6.55581e+08, mass_frac: [910.166, 2.84867e-17, 0, 0, 0]
//    PState: density: 910.166, velocity: {65.4291, -343.733, 0}, pressure: nan, temperature: nan, energy: 659071, mass_frac:
//    [1, 3.12984e-20, 0, 0, 0]
//    Previous PState: density: 910.168, velocity: {65.3781, -343.802, 0}, pressure: 9.31323e-07, temperature: 138.491, energy
//    : 659249, mass_frac: [1, 3.14697e-20, 0, 0, 0]
//
//    QState qc(910.166, {59551.3, -312854, 0}, 6.55581e+08, FractionsFlux(std::vector<double>{910.166, 0, 0, 0, 0}));
//    PState pState(qc, mixture, 9.31323e-07, 138.491);
//    std::cout << pState;
//    return;

    double max_time = 4.5e-6;
    double rho_wave = 1323.65_kg_m3, rho_water = 1000.0_kg_m3, rho_air = 1.0_kg_m3;
    double p_wave = 1.9e4_bar, p_water = 1.0_bar, p_air = 1.0_bar;
    double u_wave = 681.58_m_s, u_water = 0, u_air = 0;
    double e_wave = water->energy_rp(rho_wave, p_wave), e_water = water->energy_rp(rho_water, p_water), e_air = air->energy_rp(rho_air, p_air);
    double t_wave = water->temperature_rp(rho_wave, p_wave), t_water = water->temperature_rp(rho_water, p_water), t_air = air->temperature_rp(rho_air, p_air);

    double x_min = 0.0_cm, x_max = 0.8_cm;
    double y_min = 0.0_cm, y_max = 0.8_cm;
    double r = 0.3_cm, x_bubble = (x_max + x_min) / 2, y_bubble = (y_max + y_min) / 2;
    double x_wave = 0.06_cm;

    Fractions mass_frac_water({1, 0});
    Fractions mass_frac_air({0, 1});

    PState state_air(rho_air, Vector3d{u_air, 0, 0}, p_air, e_air, t_air, mass_frac_air);
    PState state_wave(rho_wave, Vector3d{u_wave, 0, 0}, p_wave, e_wave, t_wave, mass_frac_water);
    PState state_water(rho_water, Vector3d{u_water, 0, 0}, p_water, e_water, t_water, mass_frac_water);

    // Создаем одномерную сетку
    Rectangle rect(x_min, x_max, y_min, y_max);
    rect.set_sizes(n_cells, n_cells);
    rect.set_boundaries({
                                .left   = Boundary::ZOE, .right = Boundary::ZOE,
                                .bottom = Boundary::WALL, .top   = Boundary::WALL});

    // Файл для записи
    PvdFile pvd("mesh", filename);

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"v", get_v};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"frac1", get_frac1};
    pvd.variables += {"frac2", get_frac2};
    pvd.variables += {"c1", get_c1};
    pvd.variables += {"c2", get_c2};
    pvd.variables += {"rho1", get_rho1};
    pvd.variables += {"rho2", get_rho2};

    double time = 0.0;

    // Создать сетку
    Mesh mesh(U, &rect);

    MmFluid solver(mixture, Fluxes::GODUNOV);
    solver.set_acc(acc);

    // Число Куранта
    double CFL = 0.5;
    solver.set_CFL(CFL);

    Vector3d bubble_center = {x_bubble, y_bubble, 0};
    for (auto cell: mesh) {
        if (cell.center().x() < x_wave) {
            cell(U).set_state(state_wave);
            continue;
        }

        if ((cell.center() - bubble_center).norm() < r) {
            cell(U).set_state(state_air);
        } else {
            cell(U).set_state(state_water);
        }
    }

    double next_write = 0.0;
    int n_writes = 200;
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            std::cout << "progress: " << std::fixed << std::setprecision(1) << 100 * time / max_time << "%\n";
            pvd.save(mesh, time);
            next_write += max_time / n_writes;
        }

        // шаг решения
        solver.update(mesh);

        time = solver.get_time();
    }
}

void AirWithSF62D(int n_cells = 25, int acc = 2, const std::string &filename = "output") {
    // Уравнение состояния
    auto air = IdealGas::create("Air");
    auto sf6 = StiffenedGas::create("SF6");

    Materials mixture;
    mixture += air;
    mixture += sf6;

    double max_time = 4.0e-3;
    double rho_wave = 1153.0_kg_m3, rho_air = 1667.0_kg_m3, rho_sf6 = 5805.0_kg_m3;
    double p_wave = 1.63256e5_bar, p_air = 0.96856e5_bar, p_sf6 = 0.96856e5_bar;
    double u_wave = 133.273_m_s, u_air = 0, u_sf6 = 0;
    double e_wave = air->energy_rp(rho_wave, p_wave), e_air = air->energy_rp(rho_air, p_air), e_sf6 = sf6->energy_rp(rho_sf6, p_sf6);
    double t_wave = air->temperature_rp(rho_wave, p_wave), t_air = air->temperature_rp(rho_air, p_air), t_sf6 = sf6->temperature_rp(rho_sf6, p_sf6);

    double x_min = 0.0_m, x_max = 0.45_m;
    double y_min = 0.0_m, y_max = 0.2_m;
    double x_from = 0.1_m, x_to = 0.25_m, y_from = 0.0_m, y_to = 0.1_m;
    double x_wave = 0.02_m;

    Fractions mass_frac_air({1, 0});
    Fractions mass_frac_sf6({0, 1});

    PState state_sf6(rho_sf6, Vector3d{u_sf6, 0, 0}, p_sf6, e_sf6, t_sf6, mass_frac_sf6);
    PState state_wave(rho_wave, Vector3d{u_wave, 0, 0}, p_wave, e_wave, t_wave, mass_frac_air);
    PState state_air(rho_air, Vector3d{u_air, 0, 0}, p_air, e_air, t_air, mass_frac_air);

    // Создаем одномерную сетку
    Rectangle rect(x_min, x_max, y_min, y_max);
    rect.set_sizes(n_cells, (y_max - y_min) * n_cells / (x_max - x_min));
    rect.set_boundaries({
                                .left   = Boundary::WALL, .right = Boundary::WALL,
                                .bottom = Boundary::WALL, .top   = Boundary::WALL});

    // Файл для записи
    PvdFile pvd("mesh", filename);

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"v", get_v};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"frac1", get_frac1};
    pvd.variables += {"frac2", get_frac2};

    double time = 0.0;

    // Создать сетку
    Mesh mesh(U, &rect);

    MmFluid solver(mixture, Fluxes::GODUNOV);
    solver.set_acc(acc);

    // Число Куранта
    double CFL = 0.01;
    solver.set_CFL(CFL);

    // Настраиваем адаптацию
    mesh.set_max_level(5);
    mesh.set_distributor(solver.distributor());

    // Адаптация под начальные данные
    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        for (auto cell: mesh) {
            double x = cell.center().x(), y = cell.center().y();
            if (x < x_wave) {
                cell(U).set_state(state_wave);
                continue;
            }

            if (x_from <= x && x <= x_to && y_from <= y && y <= y_to) {
                cell(U).set_state(state_sf6);
            } else {
                cell(U).set_state(state_air);
            }
        }
        solver.set_flags(mesh);
        mesh.refine();
    }

    double next_write = 0.0;
    int n_writes = 500;
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            std::cout << "progress: " << std::fixed << std::setprecision(1) << 100 * time / max_time << "%\n";
            pvd.save(mesh, time);
            next_write += max_time / n_writes;
        }

        // шаг решения
        solver.update(mesh);

        // Установить флаги адаптации
        solver.set_flags(mesh);

        // Адаптировать сетку
        mesh.refine();

        time = solver.get_time();
    }
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

//    RiemannTesterWithSolverVertical(100, 1, "output_1");
//    RiemannTesterWithSolverVertical(100, 2, "output_3");

//    Stopwatch solve;
//    solve.start();
//    RiemannTesterWithSolver2D(Fluxes::GODUNOV, 20, 2, "output_2D_adaptive_2");
//    solve.stop();
//    std::cout << "Time: " << solve.milliseconds();
//    RiemannTesterWithSolver2D(Fluxes::GODUNOV, 20, 1, "output_2D_1");

//    KelvinHelmholtzInstability(2, "output_kelvin_3");
//    twoCellsFlux();
//    calcTests();
//    Bubble2D(30, 2, "output_bubble_4");
    Bubble2DStatic(700, 1, "output_bubble_static2");
//    AirWithSF62D(46, 2, "sf6");

    std::cout << "\nfinished\n";
    return 0;
}
