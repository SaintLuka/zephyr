#include "fast.h"
#include <zephyr/math/cfd/face_extra.h>
#include <zephyr/math/solver/mm_fluid.h>

#include <fstream>
#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/tests/classic_test.h>

#include <zephyr/math/solver/riemann.h>
#include <zephyr/phys/eos/stiffened_gas.h>
#include <filesystem>

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

double check_vacuum(AmrStorage::Item &cell) { return cell(U).p <= 0; }

struct Stat {
    double min_size, max_size, mean_size; // размеры ячеек
    int n_cells, n_min, n_max; // количество всех ячеек, количество ячеек с минимальным размером, количество ячеек с максимальным размером
    double t, dt; // прошедшее время, текущий шаг по времени
};

void write_stats_to_csv(const std::vector<Stat> &stats, const std::string &filename) {
    std::ofstream out(filename, std::ios_base::trunc);
    out << "min_size,max_size,mean_size,n_cells,n_min,n_max,time,dt\n";
    for (auto &s: stats) {
        out << s.min_size << ',' << s.max_size << ',' << s.mean_size << ','
            << s.n_cells << ',' << s.n_min << ',' << s.n_max << ','
            << s.t << ',' << s.dt
            << '\n';
    }
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

void vertical_instability_adaptive(double g = 9.81, int acc = 2, const std::string &filename = "output") {
    // Уравнение состояния
    auto mat_up = IdealGas::create(1.4, 718.0_J_kgK);
    auto mat_down = IdealGas::create(1.5, 718.0_J_kgK);

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
    double y_min = 0.0, y_max = 0.8;
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
    rect.set_sizes(8, 80);
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
        if (cell(U).mass_frac[1] == 0 && y - y_jump < 0.15 * H && abs(x - L / 2) < 0.1 * L) {
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
    std::vector<Stat> stats;
    stats.reserve(10 * n_writes);
    Stopwatch update(false), set_flags(false), refine(false);
    bool check1 = true;
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            std::cout << "progress: " << std::fixed << std::setprecision(1) << 100 * time / max_time << "%";
            std::cout << " n_cells: " << mesh.n_cells() << "\n";
            pvd.save(mesh, time);
            next_write += max_time / n_writes;
            std::cout << "Update seconds elapsed: " << update.seconds() << '\n';
            std::cout << "Set flags seconds elapsed: " << set_flags.seconds() << '\n';
            std::cout << "Refine seconds elapsed: " << refine.seconds() << "\n\n";
        }

        update.resume();
        // шаг решения
        solver.update(mesh);
        update.stop();

        double min_size = 1e6, max_size = 0, mean_size = 0;
        for (auto &cell: mesh) {
            min_size = std::min(cell.size(), min_size);
            max_size = std::max(cell.size(), max_size);
            mean_size += cell.size();
        }
        int n_min = 0, n_max = 0;
        for (auto &cell: mesh) {
            if (abs(cell.size() - min_size) / min_size < 1e-3)
                n_min++;
            else if (abs(cell.size() - max_size) / max_size < 1e-3)
                n_max++;
        }
        stats.push_back({min_size, max_size, mean_size / mesh.n_cells(), mesh.n_cells(), n_min, n_max, solver.get_time(), solver.dt()});

        set_flags.resume();
        // Установить флаги адаптации
        solver.set_flags(mesh);
        set_flags.stop();

        refine.resume();
        // Адаптировать сетку
        mesh.refine();
        refine.stop();

        if (check1 && mesh.n_cells() > 16000) {
            threads::on(16);
            check1 = false;
        }

        time = solver.get_time();
    }

    std::cout << "Update seconds elapsed: " << update.seconds() << '\n';
    std::cout << "Set flags seconds elapsed: " << set_flags.seconds() << '\n';
    std::cout << "Refine seconds elapsed: " << refine.seconds() << '\n';

    write_stats_to_csv(stats, "vertical_instability_adaptive.csv");
}

void vertical_instability_static(double g = 9.81, int acc = 2, const std::string &filename = "output") {
    // Уравнение состояния
    auto mat_up = IdealGas::create(1.4, 718.0_J_kgK);
    auto mat_down = IdealGas::create(1.5, 718.0_J_kgK);

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
    double y_min = 0.0, y_max = 0.8;
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
    rect.set_sizes(20, 500);
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

    for (auto cell: mesh) {
        if (cell.center().y() > y_jump) {
            cell(U).set_state(z_up);
        } else {
            cell(U).set_state(z_down);
        }
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

        time = solver.get_time();
    }
}

void kelvin_helmholtz_instability(int acc = 2, const std::string &filename = "output") {
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
//                cell(U).mass_frac[0] = 0.2;
//                cell(U).mass_frac[1] = 0.8;
//                cell(U).v.y() = -2;
//            }
//        }
//    }


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

void Bubble2D(int n_cells = 20, int acc = 2, const std::string &filename = "output") {
    std::filesystem::remove_all(filename); // Deletes one or more files recursively.

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

    double x_min = 0.0_cm, x_max = 1.2_cm;
    double y_min = 0.6_cm, y_max = 1.2_cm;
    double r = 0.3_cm;
    double x_bubble = 0.6_cm;
    double y_bubble = 0.6_cm;
    double x_wave = 0.06_cm;

    Fractions mass_frac_water({1, 0});
    Fractions mass_frac_air({0, 1});

    PState state_air(rho_air, Vector3d{u_air, 0, 0}, p_air, e_air, t_air, mass_frac_air);
    PState state_wave(rho_wave, Vector3d{u_wave, 0, 0}, p_wave, e_wave, t_wave, mass_frac_water);
    PState state_water(rho_water, Vector3d{u_water, 0, 0}, p_water, e_water, t_water, mass_frac_water);

    // Создаем одномерную сетку
    Rectangle rect(x_min, x_max, y_min, y_max);
    rect.set_nx(n_cells);
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
    pvd.variables += {"T", get_T};
    pvd.variables += {"frac1", get_frac1};
    pvd.variables += {"frac2", get_frac2};
    pvd.variables += {"c1", get_c1};
    pvd.variables += {"c2", get_c2};
    pvd.variables += {"rho1", get_rho1};
    pvd.variables += {"rho2", get_rho2};
    pvd.variables += {"vacuum", check_vacuum};
//    pvd.variables += {"size", [](AmrStorage::Item &cell) -> double { return cell.size; }};

    double time = 0.0;

    // Создать сетку
    Mesh mesh(U, &rect);

    MmFluid solver(mixture, Fluxes::GODUNOV);
    solver.set_acc(acc);

    // Число Куранта
    double CFL = 0.4;
    solver.set_CFL(CFL);

    // Настраиваем адаптацию
    mesh.set_max_level(6);
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
    std::vector<Stat> stats;
    stats.reserve(10 * n_writes);
    Stopwatch timer(false), set_flags(false), refine(false);
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            std::cout << "progress: " << std::fixed << std::setprecision(1) << 100 * time / max_time << "%\n";
            pvd.save(mesh, time);
            next_write += max_time / n_writes;
        }

        timer.resume();
        // шаг решения
        solver.update(mesh);
        timer.stop();

//        double min_size = 1e6, max_size = 0, mean_size = 0;
//        for (auto &cell: mesh) {
//            min_size = std::min(cell.size(), min_size);
//            max_size = std::max(cell.size(), max_size);
//            mean_size += cell.size();
//        }
//        int n_min = 0, n_max = 0;
//        for (auto &cell: mesh) {
//            if (abs(cell.size() - min_size) / min_size < 1e-3)
//                n_min++;
//            else if (abs(cell.size() - max_size) / max_size < 1e-3)
//                n_max++;
//        }
//        stats.push_back({min_size, max_size, mean_size / mesh.n_cells(), mesh.n_cells(), n_min, n_max, solver.get_time(), solver.dt()});

        set_flags.resume();
        // Установить флаги адаптации
        solver.set_flags(mesh);
        set_flags.stop();

        refine.resume();
        // Адаптировать сетку
        mesh.refine();
        refine.stop();

        time = solver.get_time();
    }

    std::cout << "Update seconds elapsed: " << timer.seconds() << '\n';
    std::cout << "Set flags seconds elapsed: " << set_flags.seconds() << '\n';
    std::cout << "Refine seconds elapsed: " << refine.seconds() << '\n';
    // Update seconds elapsed: 23
    // Set flags seconds elapsed: 5
    // Refine seconds elapsed: 108

    write_stats_to_csv(stats, "bubble2D_adaptive.csv");
}

void Bubble2DStatic(int n_cells = 100, int acc = 2, const std::string &filename = "output") {
    std::filesystem::remove_all(filename); // Deletes one or more files recursively.

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

    double x_min = 0.0_cm, x_max = 1.2_cm;
    double y_min = 0.6_cm, y_max = 1.2_cm;
    double r = 0.3_cm;
    double x_bubble = 0.6_cm;
    double y_bubble = 0.6_cm;
    double x_wave = 0.06_cm;

    Fractions mass_frac_water({1, 0});
    Fractions mass_frac_air({0, 1});

    PState state_air(rho_air, Vector3d{u_air, 0, 0}, p_air, e_air, t_air, mass_frac_air);
    PState state_wave(rho_wave, Vector3d{u_wave, 0, 0}, p_wave, e_wave, t_wave, mass_frac_water);
    PState state_water(rho_water, Vector3d{u_water, 0, 0}, p_water, e_water, t_water, mass_frac_water);

    // Создаем одномерную сетку
    Rectangle rect(x_min, x_max, y_min, y_max);
    rect.set_nx(n_cells);
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
    pvd.variables += {"T", get_T};
    pvd.variables += {"frac1", get_frac1};
    pvd.variables += {"frac2", get_frac2};
    pvd.variables += {"c1", get_c1};
    pvd.variables += {"c2", get_c2};
    pvd.variables += {"rho1", get_rho1};
    pvd.variables += {"rho2", get_rho2};
    pvd.variables += {"vacuum", check_vacuum};

    double time = 0.0;

    // Создать сетку
    Mesh mesh(U, &rect);

    MmFluid solver(mixture, Fluxes::GODUNOV);
    solver.set_acc(acc);

    // Число Куранта
    double CFL = 0.4;
    solver.set_CFL(CFL);

    Vector3d bubble_center = {x_bubble, y_bubble, 0};
    auto in_water = [bubble_center, r](const Vector3d &v) -> bool {
        return (v - bubble_center).norm() > r;
    };

    for (auto cell: mesh) {
        if (cell.center().x() < x_wave) {
            cell(U).set_state(state_wave);
            continue;
        }

        double vol_frac1 = cell.approx_vol_fraction(in_water);
        if (false && 0.0 < vol_frac1 && vol_frac1 < 1.0) {
            vol_frac1 = cell.volume_fraction(in_water, 10000);
            double vol_frac2 = 1.0 - vol_frac1;

            mmf::PState &z1 = state_water;
            mmf::PState &z2 = state_air;

            mmf::PState z(z1);

            // rho = sum a_i rho_i
            z.density = vol_frac1 * z1.density + vol_frac2 * z2.density;
            z.mass_frac = {vol_frac1 * z1.density / z.density, vol_frac2 * z2.density / z.density};
            z.velocity = z.mass_frac[0] * z1.velocity + z.mass_frac[1] * z2.velocity;
            z.pressure = z.mass_frac[0] * z1.pressure + z.mass_frac[1] * z2.pressure;

            z.temperature = mixture.temperature_rp(z.density, z.pressure, z.mass_frac);
            z.energy = mixture.energy_pt(z.pressure, z.temperature, z.mass_frac);

            cell(U).set_state(z);
        } else {
            if (vol_frac1 > 0.5) {
                cell(U).set_state(state_water);
            } else {
                cell(U).set_state(state_air);
            }
        }
    }

    double next_write = 0.0;
    int n_writes = 400;
    std::vector<Stat> stats;
    stats.reserve(10 * n_writes);
    Stopwatch timer(false);
    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            std::cout << "progress: " << std::fixed << std::setprecision(1) << 100 * time / max_time << "%\n";
            pvd.save(mesh, time);
            next_write += max_time / n_writes;
        }

        // шаг решения
        timer.resume();
        solver.update(mesh);
        timer.stop();

        double mean_size;
        for (auto &cell: mesh) {
            mean_size = cell.size();
            break;
        }
        stats.push_back({mean_size, mean_size, mean_size, mesh.n_cells(), mesh.n_cells(), mesh.n_cells(), solver.get_time(), solver.dt()});

        time = solver.get_time();
    }

    timer.stop();
    std::cout << "Update seconds elapsed: " << timer.seconds() << '\n'; // 165
//    write_stats_to_csv(stats, "bubble2D_static_test.csv");
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
    double u_wave = 10.273_m_s, u_air = 0, u_sf6 = 0;
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
    rect.set_nx(n_cells);
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
    double CFL = 0.4;
    solver.set_CFL(CFL);

    // Настраиваем адаптацию
    mesh.set_max_level(6);
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
            std::cout << "progress: " << std::fixed << std::setprecision(1) << 100 * time / max_time << "%";
            std::cout << " n_cells: " << mesh.n_cells() << "\n";
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
    threads::on(16);

//    kelvin_helmholtz_instability(2, "output_kelvin_3");
    Bubble2D(50, 2, "output_bubble_test3");
//    Bubble2DStatic(800, 1, "output_bubble_static_test");
//    AirWithSF62D(100, 2, "sf6");
//    RiemannTesterWithSolverVertical(100, 1, "output_1");
//    vertical_instability_adaptive(20, 2, "vertical_instability_adaptive2");
//    vertical_instability_static(20, 2, "vertical_instability_static");

//    Vector3d dr = Vector3d{0.00897937, 0.00620062, 0} - Vector3d{0.00898125, 0.00619875, 0};
//    std::cout << dr.x() * 99.9063 + dr.y() * (29.2343) << '\n';

    threads::off();

    return 0;
}
