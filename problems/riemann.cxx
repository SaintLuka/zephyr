#include "fast.h"

#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/tests/sod.h>
#include <zephyr/phys/tests/toro.h>

#include <zephyr/math/solver/riemann.h>
#include <zephyr/phys/eos/stiffened_gas.h>

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

using zephyr::math::RiemannSolver;

struct _U_ {
    double rho1, rho2;
    Vector3d v1, v2;
    double p1, p2;
    double e1, e2;
};

// Для быстрого доступа по типу
_U_ U;

/// Переменные для сохранения
double get_rho(Storage::Item cell) { return cell(U).rho1; }

double get_u(Storage::Item cell) { return cell(U).v1.x(); }

double get_v(Storage::Item cell) { return cell(U).v1.y(); }

double get_w(Storage::Item cell) { return cell(U).v1.z(); }

double get_p(Storage::Item cell) { return cell(U).p1; }

double get_e(Storage::Item cell) { return cell(U).e1; }


int main() {
    // Тестовая задача
    SodTest test;
//    ToroTest test(2);
    //test.inverse();

    // Уравнение состояния
    Eos& eos = test.eos;
    //StiffenedGas eos(1.367, 0.113, 0.273);

    // Состояния слева и справа в тесте
    Vector3d Ox = 100.0 * Vector3d::UnitX();
    PState zL(test.density(-Ox), test.velocity(-Ox),
              test.pressure(-Ox), test.energy(-Ox));

    PState zR(test.density(Ox), test.velocity(Ox),
              test.pressure(Ox), test.energy(Ox));

    // Точное решение задачи Римана
    RiemannSolver exact(zL, zR, eos, test.x_jump);

    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables += {"rho", get_rho};
    pvd.variables += {"u", get_u};
    pvd.variables += {"p", get_p};
    pvd.variables += {"e", get_e};

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
                      [&eos](Storage::Item cell) -> double {
                          return eos.sound_speed_rp(cell(U).rho1, cell(U).p1);
                      }};
    pvd.variables += {"c_exact",
                      [&exact, &time](const Storage::Item &cell) -> double {
                          return exact.sound_speed(cell.center().x(), time);
                      }};

    // Создаем одномерную сетку
    double H = 0.05 * (test.xmax() - test.xmin());
    Rectangle rect(test.xmin(), test.xmax(), -H, +H);
    int n_cells = 500;
    rect.set_sizes(n_cells, 1);
    rect.set_boundary_flags(
            FaceFlag::WALL, FaceFlag::WALL,
            FaceFlag::WALL, FaceFlag::WALL);

    // Создать сетку
    Mesh mesh(U, &rect);

    // Заполняем начальные данные
    for (auto cell: mesh.cells()) {
        cell(U).rho1 = test.density(cell.center());
        cell(U).v1   = test.velocity(cell.center());
        cell(U).p1   = test.pressure(cell.center());
        cell(U).e1   = eos.energy_rp(cell(U).rho1, cell(U).p1);
    }

    // Число Куранта
    double CFL = 0.9;

    // Функция вычисления потока
//    NumFlux::Ptr nf = CIR1::create();
//    NumFlux::Ptr nf = HLLC::create();
    //NumFlux::Ptr nf = Rusanov::create();
    NumFlux::Ptr nf = Godunov::create();

    double next_write = 0.0;
    size_t n_step = 0;

    while (time <= 1.01 * test.max_time()) {
        if (time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\t\tTime: " << std::setw(6) << std::setprecision(3) << time << "\n";
            pvd.save(mesh, time);
            next_write += test.max_time() / 100;
        }

        // Определяем dt
        double dt = std::numeric_limits<double>::max();
        for (auto cell: mesh.cells()) {
            // скорость звука
            double c = eos.sound_speed_rp(cell(U).rho1, cell(U).p1);
            for (auto &face: cell.faces()) {
                // Нормальная составляющая скорости
                double vn = cell(U).v1.dot(face.normal());

                // Максимальное по модулю СЗ
                double lambda = std::max(std::abs(vn + c), std::abs(vn - c));

                // Условие КФЛ
                dt = std::min(dt, cell.volume() / face.area() / lambda);
            }
        }
        dt *= CFL;

        // Расчет по некоторой схеме
        for (auto cell: mesh.cells()) {
            // Примитивный вектор в ячейке
            PState zc(cell(U).rho1, cell(U).v1, cell(U).p1, cell(U).e1);

            // Консервативный вектор в ячейке
            QState qc(zc);

            // Переменная для потока
            Flux flux;
            for (auto &face: cell.faces()) {
                // Внешняя нормаль
                auto &normal = face.normal();

                // Примитивный вектор соседа
                PState zn(zc);

                if (!face.is_boundary()) {
                    auto neib = face.neib();
                    zn.density = neib(U).rho1;
                    zn.velocity = neib(U).v1;
                    zn.pressure = neib(U).p1;
                    zn.energy = neib(U).e1;
                }

                // Значение на грани со стороны ячейки
                PState zm = zc.in_local(normal);

                // Значение на грани со стороны соседа
                PState zp = zn.in_local(normal);

                // Численный поток на грани
                auto loc_flux = nf->flux(zm, zp, eos);
                loc_flux.to_global(normal);

                // Суммируем поток
                flux.vec() += loc_flux.vec() * face.area();
            }

            // Новое значение в ячейке (консервативные переменные)
            QState Qc = qc.vec() - dt * flux.vec() / cell.volume();

            // Новое значение примитивных переменных
            PState Zc(Qc, eos);

            cell(U).rho2 = Zc.density;
            cell(U).v2   = Zc.velocity;
            cell(U).p2   = Zc.pressure;
            cell(U).e2   = Zc.energy;
        }

        // Обновляем слои
        for (auto cell: mesh.cells()) {
            std::swap(cell(U).rho1, cell(U).rho2);
            std::swap(cell(U).v1, cell(U).v2);
            std::swap(cell(U).p1, cell(U).p2);
            std::swap(cell(U).e1, cell(U).e2);
        }

        n_step += 1;
        time += dt;
    }
    
    CsvFile::save("results.csv", mesh, 5, pvd.variables);

    // расчёт ошибок
    double rho_err = 0.0, u_err = 0.0, p_err = 0.0, e_err = 0.0, c_err = 0.0;
    time = test.max_time();
    for (auto cell: mesh.cells()) {
        double x = cell.center().x();
        rho_err += abs(cell(U).rho1 - exact.density(x, time));
        u_err += abs(cell(U).v1.x() - exact.velocity(x, time));
        p_err += abs(cell(U).p1 - exact.pressure(x, time));
        e_err += abs(cell(U).e1 - exact.energy(x, time));
        c_err += abs(eos.sound_speed_rp(cell(U).rho1, cell(U).p1) - exact.sound_speed(x, time));
    }

    auto fprint = [](const std::string &name, double value) {
        std::cout << name << ": " << value << '\n';
    };
    std::cout << "Mean average errors:\n";
    fprint("\tdensity error     ", rho_err / n_cells);
    fprint("\tu error           ", u_err / n_cells);
    fprint("\tpressure error    ", p_err / n_cells);
    fprint("\tenergy error      ", e_err / n_cells);
    fprint("\tsound speed error ", c_err / n_cells);

    return 0;
}
