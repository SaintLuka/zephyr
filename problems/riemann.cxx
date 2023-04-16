#include "fast.h"

#include <zephyr/math/cfd/cir.h>
#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/tests/shu-osher.h>
#include <zephyr/phys/tests/sod.h>

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

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
    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables += {"density", get_rho};
    pvd.variables += {"velocity.x", get_u};
    pvd.variables += {"velocity.y", get_v};
    pvd.variables += {"velocity.z", get_w};
    pvd.variables += {"pressure", get_p};
    pvd.variables += {"energy", get_e};

    // Тестовая задача
    ShuOsherTest test;

    // Уравнение состояния
    Eos& eos = test.eos();

    // Создаем одномерную сетку
    double H = 0.05 * (test.xmax() - test.xmin());
    Rectangle rect(test.xmin(), test.xmax(), -H, +H);
    rect.set_sizes(1000, 1);
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
        cell(U).e1   = test.energy(cell.center());
    }

    // Число Куранта
    double CFL = 0.5;

    // Функция вычисления потока
    NumFlux::Ptr nf = CIR1::create();

    double time = 0.0;
    double next_write = 0.0;
    size_t n_step = 0;

    while (time <= 1.01 * test.max_time()) {
        if (time >= next_write) {
            std::cout << "\tШаг: " << std::setw(6) << n_step << ";"
                      << "\t\tВремя: " << std::setw(6) << std::setprecision(3) << time << "\n";
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
                    auto neib  = face.neib();
                    zn.density  = neib(U).rho1;
                    zn.velocity = neib(U).v1;
                    zn.pressure = neib(U).p1;
                    zn.energy   = neib(U).e1;
                }

                // Значение на грани со стороны ячейки
                PState zm(zc);
                zm.to_local(normal);

                // Значение на грани со стороны соседа
                PState zp(zn);
                zp.to_local(normal);

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

    return 0;
}
