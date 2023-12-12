/// @file Решение задачи распространения ударных волн в протяженных каналах с препятствиями
/// Используется схема SmFluid

#include "fast.h"
#include <zephyr/math/solver/sm_fluid.h>
#include <zephyr/phys/eos/ideal_gas.h>
#include <zephyr/math/cfd/fluxes.h>

using zephyr::math::SmFluid;
using namespace zephyr::phys;

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

EuMesh make_pipe(double xmin, double xmax, double ymin, double ymax, double H, double h, double L, double l, int nx_cells) {
    
    Rectangle rect(xmin, xmax, ymin, ymax);
    
    rect.set_nx(nx_cells);
    rect.set_boundaries({
        .left   = Boundary::WALL, .right = Boundary::WALL,
        .bottom = Boundary::WALL, .top   = Boundary::WALL});

    EuMesh mesh(U, &rect);
    
    double cell_size_x = (xmax - xmin) / nx_cells;
    double l1 = std::max(cell_size_x, l);
    
    for (auto cell: mesh) {
        if ((cell.center().y() < H-h) && (std::fmod(cell.center().x(), L) < l1)) {
            cell(U).inside = false;
        } 
        else {
            cell(U).inside = true;
        }
    }
    
    for (auto cell: mesh) {
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
    
    
    return mesh;
}


void setup_initial(EuMesh &mesh, double u0, double u3, double P0, double P3, double rho0, double rho3, double cell_size_x, IdealGas &eos) {

    for (auto cell : mesh) {
        // Инициализация
        if (cell.center().x() > 0) {
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


int main () {
    /// Начальные параметры
    /// P0, u0, rho0, P3, u3, rho3, 
    /// H, h, L, l, nx_cell (для прерятствий)

    // кПа
    double P0 = 100.0_bar;
    double P3 = 300.0_bar;

    // C
    double T0 = 20.0_C;

    // м/c
    double u0 = 0.0;
    double u3 = 299.9;

    // кг/м^3
    double rho0 = 1.1;
    double rho3 = 2.47;

    // м
    double H = 10.0; 
    double h = 5.0;
    double L = 2.0;
    double l = 0.8;

    double xmin = - 20 * L;
    double xmax = 20 * L;
    double ymin = 0.0;
    double ymax = H;

    // кол-во ячеек по x
    int nx_cells = 1000;

    // размер ячейки
    double cell_size_x = (xmax - xmin) / nx_cells;

    // сек
    double time = 0.0;
    double max_time = 1.0;
    double dt = 1.0;

    //шаг
    int n_step = 0;

    PvdFile pvd("tube", "output");

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};
    pvd.variables += {"inside", get_inside};

    // Создаем сетку
    EuMesh mesh = make_pipe(xmin, xmax, ymin, ymax, H, h, L, l, nx_cells);

    //EOS
    IdealGas eos("Air");
    Eos 

    // Инициализация начальных условий
    setup_initial(mesh, u0, u3, P0, P3, rho0, rho3, cell_size_x, eos);


    // Создать решатель
    auto solver = zephyr::math::SmFluid(mesh, flux);

    while (time <= max_time) {
        std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                  << "\tTime: " << std::setw(6) << std::setprecision(3) << time << "\n";

        pvd.save(mesh, time);

        // Шаг решения
        solver.update(mesh);

        n_step += 1;
        time += dt;
    }

    return 0;
}