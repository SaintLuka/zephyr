// Двумерное уравнение Бюргерса

#include <iostream>
#include <iomanip>

#include <zephyr/geom/boundary.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/math/funcs.h>
#include <zephyr/math/solver/burgers.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/threads.h>
#include <zephyr/utils/stopwatch.h>

using namespace zephyr::io;
using namespace zephyr::math;

using zephyr::geom::Boundary;
using zephyr::geom::generator::Rectangle;
using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::math::BurgersSolver;
using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;

// Тестовая задача
class Test {
    double angle;
    double sin_p, cos_p;
    double x_c, y_c;

public:

    Test(double alpha) : angle(alpha) {
        x_c = 0.0;
        y_c = 0.0;

        sin_p = std::sin(angle);
        cos_p = std::cos(angle);
    }

    double max_time() const { return 3.0; }

    double xmin() const { return -1.5; }
    double ymin() const { return -1.5; }
    double xmax() const { return +2.5; }
    double ymax() const { return +2.5; }

    Rectangle::Boundaries boundaries() const {
        return {.left=Boundary::ZOE, .right=Boundary::ZOE,
                .bottom=Boundary::ZOE, .top=Boundary::ZOE};
    }


    Vector3d rotate_1(const Vector3d& r) const {
        return {+cos_p * r.x() + sin_p * r.y(), -sin_p * r.x() + cos_p * r.y(), 0.0};
    }

    Vector3d rotate_2(const Vector3d& r) const {
        return {+cos_p * r.x() - sin_p * r.y(), +sin_p * r.x() + cos_p * r.y(), 0.0};
    }

    static double exact1D(double x, double t) {
        if (t == 0.0) {
            return std::abs(x) >= 1.0 ? heav(-x) : 0.5 - x + 0.5*x*std::abs(x);
        }
        if (t >= 2.0) {
            return heav(t - 2.0 * x);
        }
        if (x < t - 1.0) {
            return 1.0;
        }
        if (t - 1.0 <= x && x < 0.5 * t) {
            return (t*x + t - 1.0 + std::sqrt(2.0*t*t - 2.0*t*x - 2.0*t + 1.0)) / (t*t);
        }
        if (0.5 * t <= x && x < 1.0) {
            return (t*x - t + 1.0 - std::sqrt(2.0*t*x - 2.0*t + 1.0)) / (t*t);
        }
        return 0.0;
    }

    // Точное решение
    Vector3d operator()(const Vector3d& r, double t = 0.0) const {
        Vector3d vel = exact1D(rotate_1(r).x(), t) * Vector3d::UnitX();
        return rotate_2(vel);
    }
};

int main() {
    mpi::handler mpi_init;
    threads::on();

    Test test(0.3);

    // Генератор сетки (с граничными условиями) дает тест,
    // число ячеек можно задать
    Rectangle gen(test.xmin(), test.xmax(), test.ymin(), test.ymax());
    gen.set_boundaries(test.boundaries());
    gen.set_nx(300);

    // Создать сетку
    EuMesh mesh(gen);

    // Создать и настроить решатель
    BurgersSolver solver;
    solver.set_accuracy(2);
    solver.set_CFL(0.5);
    solver.set_limiter("MC");

    // Добавляем типы на сетку, выбираем основной слой
    auto data = solver.add_types(mesh);
    auto z = data.init;

    // Настройка сетки
    mesh.set_max_level(0);
    mesh.set_distributor(solver.distributor());

    // Задание начальных данных
    auto init_cells = [test, z](EuMesh& mesh) {
        mesh.for_each([&](EuCell& cell) {
            cell[z] = test(cell.center(), 0.0);
        });
    };

    // Файл для записи
    PvdFile pvd("Burgers", "output");

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables.append("v", z);

    double curr_time = 0.0;
    pvd.variables.append<Vector3d>("exact",
        [test, &curr_time](const EuCell &cell) -> Vector3d {
            return test(cell.center(), curr_time);
        });
    pvd.variables.append<double>("error",
        [test, z, &curr_time](const EuCell &cell) -> double {
            return (cell[z] - test(cell.center(), curr_time)).norm();
        });

    // Инициализация начальными данными
    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        init_cells(mesh);
        solver.set_flags(mesh);
        mesh.refine();
    }
    init_cells(mesh);

    size_t n_step = 0;
    double next_write = 0.0;

    Stopwatch elapsed(true);
    while (curr_time < test.max_time()) {
        if (curr_time >= next_write) {
            mpi::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(8) << std::setprecision(3) << curr_time << "\n";

            pvd.save(mesh, curr_time);
            next_write += test.max_time() / 200;
        }
        // Точное завершение в end_time
        solver.set_max_dt(test.max_time() - curr_time);

        // Обновляем слои
        solver.update(mesh);
        solver.set_flags(mesh);
        mesh.refine();

        curr_time += solver.dt();
        n_step += 1;
    }
    pvd.save(mesh, curr_time);
    elapsed.stop();

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";
    return 0;
}
