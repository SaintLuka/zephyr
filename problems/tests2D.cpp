#include "fast.h"

using zephyr::geom::generator::Rectangle;

#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/tests/sod.h>
#include <zephyr/phys/tests/toro.h>
#include <zephyr/phys/tests/blast_wave.h>
#include <zephyr/phys/tests/RiemannTest2D.h>

#include <zephyr/math/solver/riemann.h>
#include <zephyr/phys/eos/stiffened_gas.h>
#include <zephyr/math/solver/sm_fluid.h>

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

using zephyr::math::RiemannSolver;
using zephyr::math::SmFluid;

struct _U_ {
    double rho1, rho2;
    Vector3d v1, v2;
    double p1, p2;
    double e1, e2;
};

// Для быстрого доступа по типу
_U_ U;

/// Переменные для сохранения
double get_rho(AmrStorage::Item& cell) { return cell(U).rho1; }

double get_u(AmrStorage::Item& cell) { return cell(U).v1.x(); }

double get_v(AmrStorage::Item& cell) { return cell(U).v1.y(); }

double get_w(AmrStorage::Item& cell) { return cell(U).v1.z(); }

double get_p(AmrStorage::Item& cell) { return cell(U).p1; }

double get_e(AmrStorage::Item& cell) { return cell(U).e1; }


int main() {
    // Тестовая задача
    BlastWave test;
    //RiemannTest2D test(1);

    // Уравнение состояния
    Eos& eos = test.eos;
    IdealGas ig(1.4, 1.0);

    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};

    double time = 0.0;

    pvd.variables += {"c",
                      [&eos](AmrStorage::Item& cell) -> double {
                          return eos.sound_speed_rp(cell(U).rho1, cell(U).p1);
                      }};

    Rectangle gen(0, 1.0, 0, 1.0, true);
    gen.set_nx(40);
    gen.set_nx(40);
    gen.set_boundaries({.left=Boundary::WALL, .right=Boundary::WALL,
                        .bottom=Boundary::WALL, .top=Boundary::WALL});


    // Создать сетку
    EuMesh mesh(U, &gen);
    int n_cells = mesh.n_cells();

    // mesh.set_max_level(5);
        
    // Заполняем начальные данные
    Box box = mesh.bbox();
    Vector3d vc = box.center();
    for (auto cell: mesh) {
        cell(U).rho1 = test.density(cell.center() - vc);
        cell(U).v1   = test.velocity(cell.center() - vc);
        cell(U).p1   = test.pressure(cell.center() - vc);
        cell(U).e1   = eos.energy_rp(cell(U).rho1, cell(U).p1);
    }

    // Число Куранта
    double CFL = 0.5;

    // Функция вычисления потока
    // NumFlux::Ptr nf = CIR1::create();
       NumFlux::Ptr nf = HLLC::create();
    // NumFlux::Ptr nf = Rusanov::create();
    // NumFlux::Ptr nf = Godunov::create();
    // NumFlux::Ptr nf = HLL::create();

    double next_write = 0.0;
    size_t n_step = 0;

    // Создать решатель
    auto solver = zephyr::math::SmFluid();

    while (time <= test.max_time()) {
        std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                  << "\tTime: " << std::setw(6) << std::setprecision(3) << time << "\n";

        pvd.save(mesh, time);

        // Шаг решения
        solver.update(mesh, ig);

        n_step += 1;
        time += solver.m_dt;
    }

    auto fprint = [](const std::string &name, double value) {
        std::cout << name << ": " << value << '\n';
    };

    return 0;
}
