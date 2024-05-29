/// @file Решение задачи распространения ударных волн в протяженных каналах с препятствиями
/// Используется схема SmFluid

#include "fast.h"
#include <zephyr/math/solver/sm_fluid.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/math/cfd/fluxes.h>

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

struct _U_ : public SmFluid::State {
    bool inside;
};

_U_ U;

/// Переменные для сохранения
double get_rho(AmrStorage::Item& cell) { return cell(U).rho; }
double get_u(AmrStorage::Item& cell) { return cell(U).v.x(); }
double get_v(AmrStorage::Item& cell) { return cell(U).v.y(); }
double get_w(AmrStorage::Item& cell) { return cell(U).v.z(); }
double get_p(AmrStorage::Item& cell) { return cell(U).p; }
double get_e(AmrStorage::Item& cell) { return cell(U).e; }
double get_inside(AmrStorage::Item& cell) { return double(cell(U).inside); }

void setup_initial(EuMesh &mesh, double u0, double u3, double P0, double P3, double rho0, double rho3, double l, IdealGas &eos) {

    for (auto &cell: mesh) {
        // Инициализация
        if (cell.center().x() > l) {
            cell(U).v.x() = u0;
            cell(U).v.y() = 0;
            cell(U).rho = 1.0 / eos.volume_pt(P0, 20.0_C);
            cell(U).p = P0; 
            cell(U).e = eos.energy_rp(cell(U).rho, cell(U).p);
        }
        else {
            cell(U).v.x() = u3;
            cell(U).v.y() = 0;
            cell(U).rho = 1.0 / eos.volume_pt(P3, 20.0_C);
            cell(U).p = P3; 
            cell(U).e = eos.energy_rp(cell(U).rho, cell(U).p);
        }
    }
}

/// @brief 
/// @param G правая граница 
void setup_boundary(EuMesh &mesh, double G, double h, double L, double l) {

    for (auto &cell: mesh) {
        if ((cell.center().x() < G) && (cell.center().y() < h) && (std::fmod(cell.center().x(), L) < l)) {
            cell(U).inside = false;
        } 
        else {
            cell(U).inside = true;
        }
    }
    
    for (auto &cell: mesh) {
        if (!cell(U).inside) {
            continue;
        }
        for (auto &face : cell.faces()) {
            if (face.is_boundary()) {
                continue;
            }
            auto neib = face.neib();
            if (!neib(U).inside) {
                face.set_boundary(Boundary::WALL);
            }
        }
    }
}

int main () {

    std::cout << "Running with " << std::thread::hardware_concurrency() << std::endl;
    threads::on(std::thread::hardware_concurrency());

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
    switch (initials)
    {
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
        case 3:
            H = 20.0;
            h = 13.6;
            L = 4.0;
            l = 0.25;
        default:
            break;
    }

    int N = 75; // кол-во преград
    double const xmin = 0.0;
    double const xmax = N * L;
    double const ymin = 0.0;
    double const ymax = H;

    double G0 = 0.0; // положение бегущей волны 
    double G = 0.5 * L; // правая граница расчетной области
    double D0 = u3; // Скорость распространения ударной волны

    // кол-во начальных ячеек по x
    int nx_cells = 2 * N; // для 2ого случая

    // сек
    double time = 0;
    double next_write = 0;

    //шаг
    int n_step = 0;

    PvdFile pvd("tube", "/mnt/c/tube75");
    pvd.unique_nodes = true;

    //EOS
    IdealGas eos("Air");

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"inside", get_inside};
    pvd.variables += {"c",
                      [&eos](AmrStorage::Item& cell) -> double {
                          return eos.sound_speed_rp(cell(U).rho, cell(U).p);
                      }};
    pvd.variables += {"mach", 
                      [&eos](AmrStorage::Item& cell) -> double {
                          return abs(cell(U).v.x() / eos.sound_speed_rp(cell(U).rho, cell(U).p));
                      }};

    // Создаем сетку
    Rectangle rect(xmin, xmax, ymin, ymax);
    rect.set_nx(nx_cells);
    rect.set_boundaries({
        .left   = Boundary::ZOE, .right = Boundary::ZOE,
        .bottom = Boundary::WALL, .top   = Boundary::WALL});

    EuMesh mesh(U, &rect);

    // Создать решатель
    auto solver = zephyr::math::SmFluid(eos, Fluxes::HLLC);
    solver.set_accuracy(2);
    solver.set_CFL(0.4);

    mesh.set_max_level(3);
    mesh.set_distributor(solver.distributor());

    for (int k = 0; k < mesh.max_level() + 1; ++k) {

        setup_initial(mesh, u0, u3, P0, P3, rho0, rho3, l, eos);

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
    setup_boundary(mesh, G, h, L, l);

    pvd.save(mesh, time);

    while (1) {

        std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                  << "\tTime: " << std::setw(6) << std::setprecision(3) << time << "\n";

        if (time >= next_write) {
            pvd.save(mesh, time);
            next_write += 0.01;
        };

        // Обновляем слои
        solver.update(mesh);

        // TODO
        G0 += 1.4 * D0 * solver.get_m_dt();

        if (G - G0 < L) {

            G += L;

            // обновить сетку
            for (int k = 0; k < mesh.max_level() + 1; ++k) {
                for (auto &cell: mesh) {
                    bool need_split = false;
                    if (cell.center().x() < G) {
                        need_split = true;    
                    }
                    if(need_split) {
                        cell.set_flag(1);
                    } else {
                        cell.set_flag(-1);
                    }
                }
                mesh.refine(); 
            }

            setup_boundary(mesh, G, h, L, l);
        }
        
        n_step += 1;
        time = solver.get_time();
    }

    auto fprint = [](const std::string &name, double value) {
        std::cout << name << ": " << value << '\n';
    };

    threads::off();

    return 0;
}