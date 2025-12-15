/// @file clogged_shock_tube.cpp
/// @brief Распространение ударных волн в протяженных каналах с препятствиями.
/// Используется решатель SmFluid.
#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/utils/threads.h>
#include <zephyr/io/pvd_file.h>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/solver/sm_fluid.h>

using namespace zephyr::mesh;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;
using namespace zephyr::io;

using generator::Rectangle;
using zephyr::utils::threads;

void setup_initial(EuMesh &mesh, Storable<PState> z, IdealGas &eos, double u0, double u3, double P0, double P3, double rho0, double rho3, double l) {
    for (auto cell: mesh) {
        // Инициализация
        if (cell.center().x() > l) {
            cell[z].velocity.x() = u0;
            cell[z].velocity.y() = 0;
            cell[z].density = 1.0 / eos.volume_PT(P0, 20.0_C);
            cell[z].pressure = P0;
        }
        else {
            cell[z].velocity.x() = u3;
            cell[z].velocity.y() = 0;
            cell[z].density = 1.0 / eos.volume_PT(P3, 20.0_C);
            cell[z].pressure = P3;
        }
        cell[z].energy = eos.energy_rP(cell[z].density, cell[z].pressure);
    }
}

/// @param G Правая граница
void setup_boundary(EuMesh &mesh, Storable<bool> inside, double G, double h, double L, double l) {
    for (auto &cell: mesh) {
        cell[inside] = !((cell.center().x() < G) &&
                         (cell.center().y() < h) &&
                         (std::fmod(cell.center().x(), L) < l));
    }
    
    for (auto cell: mesh) {
        if (!cell[inside]) {
            continue;
        }
        for (auto face: cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }
            if (!face.neib(inside)) {
                face.set_boundary(Boundary::WALL);
            }
        }
    }
}

int main() {
    threads::on();

    int initials = 3;
    // 0 - 4.2.1 Постановка задачи                                      Таблица 3
    // 1 - 4.2.2 Эффективная ударная адиабата канала с препятствиями    Таблица 4
    // 2 - 4.2.4 Течение газа в сильно загромождённых каналах           Таблица 5

    /// Начальные параметры
    /// P0, u0, rho0, P3, u3, rho3, 
    /// H, h, L, l, nx_cell (для препятствий)

    // кПа
    double P0 = 100.0_bar;
    double P3 = 300.0_bar;

    // C
    double T0 = 20.0_C;

    // м/c
    double u0 = 0.0;
    double u3 = 299.9;

    // кг/м^3
    double rho0 = 1.17;
    double rho3 = 2.47;

    // м
    double H, h, L, l;
    switch (initials) {
        case 0:
            H = 10.0; 
            h = 5.0;
            L = 2.0;
            l = 0.8;
            break;
        case 1:
            H = 10.0; 
            h = 5.0;
            L = 4.0;
            l = 0.8;
            break;
        case 2:
            H = 20.0;
            h = 13.6;
            L = 4.0;
            l = 0.2;
            break;
        case 3:
            H = 20.0;
            h = 13.6;
            L = 4.0;
            l = 0.25;
            break;
        default:
            throw std::runtime_error("Wrong initials");
    }

    int N = 75; // кол-во преград
    double const xmin = 0.0;
    double const xmax = N * L;
    double const ymin = 0.0;
    double const ymax = H;

    double G0 = 0.0;    // положение бегущей волны
    double G = 0.5 * L; // правая граница расчетной области
    double D0 = u3;     // Скорость распространения ударной волны

    // кол-во начальных ячеек по x
    int nx_cells = 2 * N; // для 2ого случая

    // сек
    double curr_time = 0;
    double next_write = 0;

    //шаг
    int n_step = 0;

    // EOS
    IdealGas::Ptr eos = IdealGas::create("Air");

    // Генератор сетки
    Rectangle rect(xmin, xmax, ymin, ymax);
    rect.set_nx(nx_cells);
    rect.set_boundaries({
        .left   = Boundary::ZOE,  .right = Boundary::ZOE,
        .bottom = Boundary::WALL, .top   = Boundary::WALL});

    // Создать сетку
    EuMesh mesh(rect);

    // Создать решатель
    SmFluid solver(eos);
    solver.set_accuracy(2);
    solver.set_CFL(0.4);
    solver.set_method(Fluxes::HLLC);

    // Добавляем типы на сетку, выбираем основной слой
    auto data = solver.add_types(mesh);
    auto z = data.init;
    auto inside = mesh.add<bool>("inside");

    // Настройка сетки
    mesh.set_max_level(3);
    mesh.set_distributor(solver.distributor());

    PvdFile pvd("tube", "output");
    //pvd.unique_nodes = true;

    // Переменные для сохранения
    pvd.variables = {"level", "faces2D"};
    pvd.variables += {"rho", [z](EuCell& cell) -> double { return cell[z].density; }};
    pvd.variables += {"u",   [z](EuCell& cell) -> double { return cell[z].velocity.x(); }};
    pvd.variables += {"v",   [z](EuCell& cell) -> double { return cell[z].velocity.y(); }};
    pvd.variables += {"p",   [z](EuCell& cell) -> double { return cell[z].pressure; }};
    pvd.variables += {"e",   [z](EuCell& cell) -> double { return cell[z].energy; }};
    pvd.variables += {"inside", [inside](EuCell& cell) -> int { return cell[inside]; } };
    pvd.variables += {"c", [&eos, z](EuCell& cell) -> double {
        return eos->sound_speed_rP(cell[z].density, cell[z].pressure);
    }};
    pvd.variables += {"mach", [&eos, z](EuCell& cell) -> double {
        return abs(cell[z].velocity.x() / eos->sound_speed_rP(cell[z].density, cell[z].pressure));
    }};

    for (int k = 0; k < mesh.max_level() + 1; ++k) {

        setup_initial(mesh, z, *eos, u0, u3, P0, P3, rho0, rho3, l);

        for (auto &cell: mesh) {
            bool need_split = false;
            if (cell.center().x() < G) {
                need_split = true;
            }
            if (need_split) {
                cell.set_flag(1);
            } else {
                cell.set_flag(-1);
            }
        }

        mesh.refine();
    }

    // Инициализация граничных условий
    setup_boundary(mesh, inside, G, h, L, l);

    while (curr_time < 0.8) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";

            pvd.save(mesh, curr_time);
            next_write += 0.01;
        }

        // Обновляем слои
        solver.update(mesh);

        // TODO
        G0 += 1.4 * D0 * solver.dt();

        if (G - G0 < L) {
            G += L;

            // обновить сетку
            for (int k = 0; k < mesh.max_level() + 1; ++k) {
                for (auto &cell: mesh) {
                    bool need_split = false;
                    if (cell.center().x() < G) {
                        need_split = true;    
                    }
                    if (need_split) {
                        cell.set_flag(1);
                    } else {
                        cell.set_flag(-1);
                    }
                }
                mesh.refine(); 
            }

            setup_boundary(mesh, inside, G, h, L, l);
        }

        curr_time += solver.dt();
        n_step += 1;
    }

    return 0;
}