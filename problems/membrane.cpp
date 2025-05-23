/// @file membrane.cpp
/// @brief Решатель газодинамики в одном файле.
/// В углу прямоугольной области колеблется мембрана.

#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/rectangle.h>

#include <zephyr/mesh/euler/soa_mesh.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/math/cfd/fluxes.h>
#include <zephyr/math/cfd/models.h>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/matter/eos/ideal_gas.h>

#include <zephyr/utils/stopwatch.h>

using zephyr::geom::Box;
using zephyr::geom::Boundary;
using zephyr::geom::Vector3d;
using zephyr::mesh::generator::Rectangle;
using zephyr::mesh::SoaMesh;
using zephyr::mesh::QCell;
using zephyr::io::PvdFile;
using zephyr::utils::Stopwatch;
using zephyr::utils::threads;

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::smf;

int main() {
    // Включаем многопоточность
    threads::on();

    // Тестовая задача
    IdealGas eos("Air");
    double T0 = 20.0_C;
    double P0 = 1.0_bar;
    double R0 = 1.0 / eos.volume_PT(P0, T0);

    // Создаем одномерную сетку
    Rectangle rect(0.0, 4.0, -2.0, 4.0, true);
    rect.set_nx(200);
    rect.set_boundaries({
                                .left   = Boundary::WALL, .right = Boundary::WALL,
                                .bottom = Boundary::WALL, .top   = Boundary::WALL});

    // Создать сетку
    SoaMesh mesh(rect);

    // Переменные для хранения на сетке
    auto [rho1, p1, e1] = mesh.append<double>("rho1", "p1", "e1");
    auto [rho2, p2, e2] = mesh.append<double>("rho2", "p2", "e2");
    auto [v1, v2]       = mesh.append<Vector3d>("v1", "v2");

    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables += {"rho", [rho1](QCell& cell) -> double { return cell(rho1); }};
    pvd.variables += {"velocity.x", [v1](QCell& cell) -> double { return cell(v1).x(); }};
    pvd.variables += {"velocity.y", [v1](QCell& cell) -> double { return cell(v1).y(); }};
    pvd.variables += {"|velocity|", [v1](QCell& cell) -> double { return cell(v1).norm(); }};
    pvd.variables += {"pressure",   [p1](QCell& cell) -> double { return cell(p1); }};
    pvd.variables += {"SPL", [p1](QCell& cell) -> double {
        // Уровень звукового давления в дБ
        return 20.0 * std::log(1.0 + std::abs(cell(p1) - 1.0_bar) / 20.0e-6_Pa) / std::log(10.0);
    }};

    // Заполняем начальные данные
    for (auto cell: mesh) {
        cell(rho1) = R0;
        cell(v1) = Vector3d::Zero();
        cell(p1) = P0;
        cell(e1) = eos.energy_rP(cell(rho1), cell(p1));
    }

    // Число Куранта
    double CFL = 0.5;

    //NumFlux::Ptr nf = CIR1::create();
    NumFlux::Ptr nf = HLLC::create();

    double time = 0.0;
    double next_write = 0.0;
    size_t n_step = 0;
    double max_time = 0.05;
    double write_freq = max_time / 100;

    Stopwatch elapsed(true);
    Stopwatch write;
    Stopwatch sw_dt;
    Stopwatch sw_flux;
    Stopwatch sw_update;

    while (time <= 1.01 * max_time) {
        if (time >= next_write) {
            write.resume();
            std::cout << "\tStep: " << std::setw(8) << n_step << ";\t\t"
                      << "Time: " << std::setw(11) << std::setprecision(3) << time << "\n";
            pvd.save(mesh, time);
            next_write += write_freq;
            write.stop();
        }

        // Определяем dt
        sw_dt.resume();
        double dt = mesh.min([&](QCell cell) -> double {
            double dt = 1.0e300;

            // скорость звука
            double c = eos.sound_speed_rP(cell(rho1), cell(p1));
            for (auto &face: cell.faces()) {
                // Нормальная составляющая скорости
                double vn = cell(v1).dot(face.normal());

                // Максимальное по модулю СЗ
                double lambda = std::max(std::abs(vn + c), std::abs(vn - c));

                // Условие КФЛ
                dt = std::min(dt, cell.volume() / face.area() / lambda);
            }

            return dt;
        }, 1.0e300);
        dt *= CFL;
        sw_dt.stop();

        // Расчет по схеме CIR
        sw_flux.resume();
        mesh.for_each([&](QCell cell) {
            // Примитивный вектор в ячейке
            PState zc(cell(rho1), cell(v1), cell(p1), cell(e1));

            // Консервативный вектор в ячейке
            QState qc(zc);

            // Переменная для потока
            Flux flux;
            for (auto &face: cell.faces()) {
                // Внешняя нормаль
                auto& normal = face.normal();

                // Примитивный вектор соседа
                PState zn(zc);

                if (face.is_boundary()) {
                    const double d = 0.5;

                    // Граничные условия типа "стенка"
                    double vn = zc.velocity.dot(face.normal());
                    zn.velocity -= 2.0 * vn * face.normal();

                    // Камертон
                    if (face.normal().x() < 0.0 && fabs(face.y()) < d) {
                        // Амплитуда колебаний
                        double A = 0.0;
                        if (fabs(face.y()) < d) {
                            A = std::pow(std::cos(M_PI_2 * face.y() / d), 2);
                        }

                        // Непосредственно колебание
                        double T = std::sin(2.0 * M_PI * 440.0 * time);

                        // Затухание
                        double R = std::exp(-100.0 * time);

                        zn.velocity.x() += 1.0e-5 * R * T * A;
                    }
                } else {
                    zn.density  = face.neib(rho1);
                    zn.velocity = face.neib(v1);
                    zn.pressure = face.neib(p1);
                    zn.energy   = face.neib(e1);
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

            cell(rho2) = Zc.density;
            cell(v2)   = Zc.velocity;
            cell(p2)   = Zc.pressure;
            cell(e2)   = Zc.energy;
        });
        sw_flux.stop();

        // Обновляем слои
        sw_update.resume();
        mesh.swap(rho1, rho2);
        mesh.swap(v1, v2);
        mesh.swap(p1, p2);
        mesh.swap(e1, e2);
        sw_update.stop();

        n_step += 1;
        time += dt;
    }
    elapsed.stop();

    std::cout << "\nElapsed time: " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    std::cout << "  Write time:   " << write.extended_time()
              << " ( " << write.milliseconds() << " ms)\n";

    std::cout << "  dt time:      " << sw_dt.extended_time()
              << " ( " << sw_dt.milliseconds() << " ms)\n";

    std::cout << "  flux time:    " << sw_flux.extended_time()
              << " ( " << sw_flux.milliseconds() << " ms)\n";

    std::cout << "  update time:  " << sw_update.extended_time()
              << " ( " << sw_update.milliseconds() << " ms)\n";

    return 0;
}