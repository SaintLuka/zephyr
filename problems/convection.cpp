/// @file convection.cpp
/// @brief Решение задачи переноса с неоднородной скоростью.
/// Используется схема высокого порядка (MUSCL) с адаптацией.
#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/cuboid.h>

#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/math/solver/convection.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/utils/threads.h>
#include <zephyr/utils/stopwatch.h>

using zephyr::geom::Box;
using zephyr::geom::Boundary;
using zephyr::geom::Vector3d;
using zephyr::geom::generator::Rectangle;
using zephyr::geom::generator::Cuboid;
using zephyr::mesh::EuCell;
using zephyr::mesh::EuMesh;
using zephyr::mesh::Storable;
using zephyr::io::PvdFile;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;


/// @brief Наследуем собственный решатель от Convection, теперь переопределив
/// поле скорости можно решать произвольные задачи на перенос.
class Solver : public zephyr::math::Convection {
public:

    /// @brief Задаем скорость переноса
    Vector3d velocity(const Vector3d& c) const override {
        //return { 0.5, 0.4, 0.3 };
        return { 0.8, 0.40 + 0.3*std::sin(4 * M_PI * c.x()), 0.2 };
    }
};

const double margin = 0.0198;

// Начальное условие в виде круга
void setup_initial_1(EuMesh& mesh, double D, Storable<double> u) {
    //double R = D / 2.0;
    //Vector3d vc = {R + margin, R + margin, 0.2};
    double R = 0.1;
    Vector3d vc = {0.15, 0.15, 0.0};
    for (auto cell: mesh) {
        cell(u) = (cell.center() - vc).norm() < R ? 1.0 : 0.0;
    }
}

// Начальное условие в виде квадрата
void setup_initial_2(EuMesh& mesh, double D, Storable<double> u) {
    double x_min = margin;
    double x_max = D + x_min;
    double y_min = margin;
    double y_max = D + y_min;

    for (auto cell: mesh) {
        Vector3d vc = cell.center();
        if (x_min <= vc.x() && vc.x() <= x_max &&
            y_min <= vc.y() && vc.y() <= y_max) {
            cell(u) = 1.0;
        } else {
            cell(u) = 0.0;
        }
        cell(u) = 0.0;
    }
}

int main() {
    threads::on();

    // Геометрия области
#define GENTYPE 3

#if GENTYPE == 1
    // Прямоугольник
    Rectangle gen(0.0, 1.0, 0.0, 0.6, false);
    gen.set_nx(100);
    gen.set_boundaries({
        .left   = Boundary::ZOE, .right = Boundary::ZOE,
        .bottom = Boundary::ZOE, .top   = Boundary::ZOE});

#elif GENTYPE == 2
    // Прямоугольник из сот
    Rectangle gen(0.0, 1.0, 0.0, 0.6, true);
    gen.set_nx(500);
    gen.set_boundaries({
        .left   = Boundary::ZOE, .right = Boundary::ZOE,
        .bottom = Boundary::ZOE, .top   = Boundary::ZOE});

#elif GENTYPE == 3
    // Кубоид
    Cuboid gen(0.0, 1.0, 0.0, 0.8, -0.3, 0.3);
    gen.set_nx(20);
    gen.set_boundaries({
        .left   = Boundary::ZOE, .right = Boundary::ZOE,
        .bottom = Boundary::ZOE, .top   = Boundary::ZOE,
        .back   = Boundary::ZOE, .front = Boundary::ZOE});
#endif

    // Создать решатель
    Solver solver;

    // Настроим решатель
    solver.set_CFL(0.5);
    solver.set_accuracy(3);
    solver.set_limiter("MC");

    // Создать сетку
    EuMesh mesh(gen);

    /*
    if (mesh.check_base() < 0) {
        int res = mesh.check_base();
        std::cout << "Bad init mesh " << res << "\n";
        return 0;
    }
     */

    // Добавить переменные на сетку
    solver.add_types(mesh);

    // Настраиваем адаптацию
    mesh.set_max_level(mesh.dim() < 3 ? 5 : 4);
    mesh.set_distributor(solver.distributor());

    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables += {"u", [&solver](EuCell& cell) { return cell(solver.u_curr); } };
    pvd.variables += {"u2", [&solver](EuCell& cell) { return cell(solver.u_next); } };
    pvd.variables += {"uh", [&solver](EuCell& cell) { return cell(solver.u_half); } };
    pvd.variables += {"dx", [&solver](EuCell& cell) { return cell(solver.du_dx); } };
    pvd.variables += {"dy", [&solver](EuCell& cell) { return cell(solver.du_dy); } };
    pvd.variables += {"vx", [&solver](EuCell& cell) { return solver.velocity(cell.center()).x(); } };
    pvd.variables += {"vy", [&solver](EuCell& cell) { return solver.velocity(cell.center()).y(); } };

    // Заполняем начальные данные
    Box box = mesh.bbox();
    double D = 0.1*box.diameter();

    // Адаптация под начальные данные
    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        setup_initial_1(mesh, D, solver.u_curr);
        solver.set_flags(mesh);
        mesh.refine();

        /*
        if (mesh.check_refined() < 0) {
            int res = mesh.check_refined();
            std::cout << "Bad init refinement\n";
            return 0;
        }
         */
    }

    Stopwatch elapsed;
    Stopwatch sw_write;
    Stopwatch sw_update;
    Stopwatch sw_flags;
    Stopwatch sw_refine;

    int n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    elapsed.start();
    while(curr_time <= 1.0) {
        sw_write.resume();
        if (curr_time >= next_write) {
            std::cout << "\tШаг: " << std::setw(6) << n_step << ";"
                      << "\tВремя: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);
            next_write += 0.02;
        }
        sw_write.stop();

        // Шаг решения
        sw_update.resume();
        solver.update(mesh);
        sw_update.stop();

        // Установить флаги адаптации
        sw_flags.resume();
        solver.set_flags(mesh);
        sw_flags.stop();

        // Адаптировать сетку
        sw_refine.resume();
        mesh.refine();
        sw_refine.stop();

        n_step += 1;
        curr_time += solver.dt();
    }
    elapsed.stop();

    std::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    std::cout << "  Write time:   " << sw_write.extended_time()
              << " ( " << sw_write.milliseconds() << " ms)\n";

    std::cout << "  Update time:  " << sw_update.extended_time()
              << " ( " << sw_update.milliseconds() << " ms)\n";

    std::cout << "  Set flags:    " << sw_flags.extended_time()
              << " ( " << sw_flags.milliseconds() << " ms)\n";

    std::cout << "  Refinement:   " << sw_refine.extended_time()
              << " ( " << sw_refine.milliseconds() << " ms)\n";

    return 0;
}