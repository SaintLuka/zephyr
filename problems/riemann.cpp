#include "fast.h"

#include <zephyr/geom/generator/collection/wedge.h>
#include <zephyr/geom/generator/collection/semicircle_cutout.h>

using zephyr::geom::generator::collection::Wedge;
using zephyr::geom::generator::collection::SemicircleCutout;

#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/tests/sod.h>
#include <zephyr/phys/tests/toro.h>
#include <zephyr/phys/tests/shu_osher.h>

#include <zephyr/math/solver/riemann.h>
#include <zephyr/phys/eos/stiffened_gas.h>
#include <zephyr/math/solver/sm_fluid.h>

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

using zephyr::math::RiemannSolver;


// Для быстрого доступа по типу
SmFluid::State U;

/// Переменные для сохранения
double get_rho(AmrStorage::Item& cell) { return cell(U).rho; }
double get_u(AmrStorage::Item& cell) { return cell(U).v.x(); }
double get_v(AmrStorage::Item& cell) { return cell(U).v.y(); }
double get_w(AmrStorage::Item& cell) { return cell(U).v.z(); }
double get_p(AmrStorage::Item& cell) { return cell(U).p; }
double get_e(AmrStorage::Item& cell) { return cell(U).e; }


int main() {
    // Тестовая задача
    // SodTest test;
    // ToroTest test(3);
    // test.inverse();
    ShuOsherTest test;

    // Уравнение состояния
    Eos& eos = test.eos();
    //StiffenedGas eos(1.367, 0.113, 0.273);
    StiffenedGas sg = eos.stiffened_gas(1.0, 1.0);

    // Состояния слева и справа в тесте
    Vector3d Ox = 100.0 * Vector3d::UnitX();
    PState zL(test.density(-Ox), test.velocity(-Ox),
              test.pressure(-Ox), test.energy(-Ox));

    PState zR(test.density(Ox), test.velocity(Ox),
              test.pressure(Ox), test.energy(Ox));

    // Точное решение задачи Римана
    //RiemannSolver exact(zL, zR, sg, test.x_jump);
    // RiemannSolver exact(zL, zR, sg, -4.0);

    // Файл для записи
    PvdFile pvd("mesh", "/mnt/d/ShuOsher400"); 

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};

    double time = 0.0;

    // pvd.variables += {"rho_exact",
    //                   [&exact, &time](const AmrStorage::Item &cell) -> double {
    //                       return exact.density(cell.center.x(), time);
    //                   }};
    // pvd.variables += {"u_exact",
    //                   [&exact, &time](const AmrStorage::Item &cell) -> double {
    //                       return exact.velocity(cell.center.x(), time);
    //                   }};
    // pvd.variables += {"p_exact",
    //                   [&exact, &time](const AmrStorage::Item &cell) -> double {
    //                       return exact.pressure(cell.center.x(), time);
    //                   }};
    // pvd.variables += {"e_exact",
    //                   [&exact, &time](const AmrStorage::Item &cell) -> double {
    //                       return exact.energy(cell.center.x(), time);
    //                   }};
    // pvd.variables += {"c",
    //                   [&eos](AmrStorage::Item& cell) -> double {
    //                       return eos.sound_speed_rp(cell(U).rho, cell(U).p);
    //                   }};
    // pvd.variables += {"c_exact",
    //                   [&exact, &time](const AmrStorage::Item &cell) -> double {
    //                       return exact.sound_speed(cell.center.x(), time);
    //                   }};

    // Часть области с регулярной сеткой
    // auto fix_condition = [&test](const Vector3d& v) {
    //     return v.x() <= test.x_jump + 0.01 * (test.xmax() - test.xmin());
    // };

    // Создаем одномерную сетку
    Strip gen(test.xmin(), test.xmax());
    int n_cells = 400;
    gen.set_size(n_cells);

    //Wedge gen(0.40, 0.9, 0.0, 0.20, 0.6, 0.1 * M_PI);
    //gen.set_nx(200);
    //gen.set_fixed(fix_condition);
    //gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::WALL,
    //                    .bottom=Boundary::WALL, .top=Boundary::ZOE});

    // SemicircleCutout gen(0.49, 0.7, 0.0, 0.07, 0.6, 0.02);
    // gen.set_ny(10);
    // gen.set_fixed(fix_condition);
    // gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
    //                     .bottom=Boundary::WALL, .top=Boundary::ZOE});

    // Создать сетку
    EuMesh mesh(U, &gen);
    //int n_cells = mesh.n_cells();

    // Создать решатель
    auto solver = zephyr::math::SmFluid(eos, Fluxes::HLLC);

    //mesh.set_max_level(5);
    //mesh.set_distributor(solver.distributor());
    solver.init_cells(mesh, test);
    solver.set_accuracy(2);

    // Число Куранта
    double CFL = 0.2;
    solver.set_CFL(CFL);

    double next_write = 0.0;
    size_t n_step = 0;

    while (time <= 5000.01 * test.max_time()) {
        if (time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << time << "\n";
            pvd.save(mesh, time);
            next_write += test.max_time() * 20;
        }

        // Обновляем слои
        solver.update(mesh);

        n_step += 1;
        time += solver.get_time();
    }

    auto fprint = [](const std::string &name, double value) {
        std::cout << name << ": " << value << '\n';
    };

    return 0;
}
