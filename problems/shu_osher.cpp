/// @file shu_oscher.cpp
/// @brief Решатель газодинамики в одном файле на примере теста Шу-Ошера

#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/strip.h>

#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>

#include <zephyr/phys/tests/test_1D.h>

using zephyr::geom::Box;
using zephyr::geom::Boundary;
using zephyr::geom::Vector3d;
using zephyr::geom::generator::Strip;
using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::io::PvdFile;

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;


int main() {
    // Тестовая задача
    ShuOsherTest test;

    // Уравнение состояния
    auto eos = test.get_eos();

    // Создаем одномерную сетку
    Strip gen(test.xmin(), test.xmax());
    gen.set_size(10000);
    gen.set_boundaries({.left = Boundary::ZOE, .right = Boundary::ZOE});

    // Создать сетку
    EuMesh mesh(gen);

    // Переменные для хранения на сетке
    auto [rho1, p1, e1] = mesh.add<double>("rho1", "p1", "e1");
    auto [rho2, p2, e2] = mesh.add<double>("rho2", "p2", "e2");
    auto [v1, v2]       = mesh.add<Vector3d>("v1", "v2");

    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables += {"rho",    [rho1](EuCell& cell) -> double { return cell(rho1); }};
    pvd.variables += {"velocity", [v1](EuCell& cell) -> double { return cell(v1).x(); }};
    pvd.variables += {"pressure", [p1](EuCell& cell) -> double { return cell(p1); }};
    pvd.variables += {"energy",   [e1](EuCell& cell) -> double { return cell(e1); }};

    // Заполняем начальные данные
    for (auto cell: mesh) {
        cell(rho1) = test.density (cell.center());
        cell(v1)   = test.velocity(cell.center());
        cell(p1)   = test.pressure(cell.center());
        cell(e1)   = test.energy  (cell.center());
    }

    // Число Куранта
    double CFL = 0.5;

    // Функция вычисления потока
    NumFlux::Ptr nf = HLLC::create();

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
        for (auto cell: mesh) {
            // скорость звука
            double c = eos->sound_speed_rP(cell(rho1), cell(p1));
            for (auto &face: cell.faces()) {
                // Нормальная составляющая скорости
                double vn = cell(v1).dot(face.normal());

                // Максимальное по модулю СЗ
                double lambda = std::max(std::abs(vn + c), std::abs(vn - c));

                // Условие КФЛ
                dt = std::min(dt, cell.volume() / face.area() / lambda);
            }
        }
        dt *= CFL;

        // Расчет по некоторой схеме
        for (auto cell: mesh) {
            // Примитивный вектор в ячейке
            PState zc(cell(rho1), cell(v1), cell(p1), cell(e1));

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
                    zn.density  = face.neib(rho1);
                    zn.velocity = face.neib(v1);
                    zn.pressure = face.neib(p1);
                    zn.energy   = face.neib(e1);
                }

                // Значение на грани со стороны ячейки
                PState zm(zc);
                zm.to_local(normal);

                // Значение на грани со стороны соседа
                PState zp(zn);
                zp.to_local(normal);

                // Численный поток на грани
                auto loc_flux = nf->flux(zm, zp, *eos);
                loc_flux.to_global(normal);

                // Суммируем поток
                flux.vec() += loc_flux.vec() * face.area();
            }

            // Новое значение в ячейке (консервативные переменные)
            QState Qc = qc.vec() - dt * flux.vec() / cell.volume();

            // Новое значение примитивных переменных
            PState Zc(Qc, *eos);

            cell(rho2) = Zc.density;
            cell(v2)   = Zc.velocity;
            cell(p2)   = Zc.pressure;
            cell(e2)   = Zc.energy;
        }

        // Обновляем слои
        mesh.swap(rho1, rho2);
        mesh.swap(v1, v2);
        mesh.swap(p1, p2);
        mesh.swap(e1, e2);

        n_step += 1;
        time += dt;
    }

    return 0;
}