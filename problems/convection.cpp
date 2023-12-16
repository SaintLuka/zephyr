/// @file Решение задачи переноса с неоднородной скоростью.
/// Используется схема высокого порядка (MUSCL) с адаптацией.

#include "fast.h"

#include <zephyr/math/solver/convection.h>

/// @brief Наследуем собственный решатель от Convection, теперь переопределив
/// поле скорости можно решать произвольные задачи на перенос.
class Solver : public zephyr::math::Convection {
public:

    /// @brief Задаем скорость переноса
    Vector3d velocity(const Vector3d& c) const override {
        return { 0.5, 0.4, 0.3 };
        return { 0.8, 0.40 + 0.3*std::sin(4 * M_PI * c.x()), 0.2 };
    }
};

// Получим тип данных
Solver::State U = Solver::datatype();

// Переменные для сохранения
double get_u(AmrStorage::Item& cell)  { return cell(U).u1; }

double get_ux(AmrStorage::Item& cell)  { return cell(U).ux; }

double get_uy(AmrStorage::Item& cell)  { return cell(U).uy; }

double get_lvl(AmrStorage::Item& cell)  { return cell.level; }

const double margin = 0.0198;

// Начальное условие в виде круга
void setup_initial_1(Mesh& mesh, double D) {
    //double R = D / 2.0;
    //Vector3d vc = {R + margin, R + margin, 0.2};
    double R = 0.1;
    Vector3d vc = {0.15, 0.15, 0.0};
    for (auto cell: mesh) {
        cell(U).u1 = (cell.center() - vc).norm() < R ? 1.0 : 0.0;
        cell(U).u2 = 0.0;
    }
}

// Начальное условие в виде квадрата
void setup_initial_2(Mesh& mesh, double D) {
    double x_min = margin;
    double x_max = D + x_min;
    double y_min = margin;
    double y_max = D + y_min;

    for (auto cell: mesh) {
        Vector3d vc = cell.center();
        if (x_min <= vc.x() && vc.x() <= x_max &&
            y_min <= vc.y() && vc.y() <= y_max) {
            cell(U).u1 = 1.0;
        } else {
            cell(U).u1 = 0.0;
        }
        cell(U).u2 = 0.0;
    }
}

int main() {
    threads::on();

    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables += {"u",  get_u};
    pvd.variables += {"ux", get_ux};
    pvd.variables += {"uy", get_uy};
    pvd.variables += {"lvl", get_lvl};

    // Геометрия области
#define GENTYPE 1

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
    gen.set_nx(300);
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
    solver.set_limiter("van Leer");

    // Создать сетку
    Mesh mesh(U, &gen);

    // Настраиваем адаптацию
    mesh.set_max_level(5);
    mesh.set_distributor(solver.distributor());

    // Заполняем начальные данные
    Box box = mesh.bbox();
    double D = 0.1*box.diameter();

    // Адаптация под начальные данные
    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        setup_initial_1(mesh, D);
        solver.set_flags(mesh);
        mesh.refine();
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
            pvd.save(mesh.locals(), curr_time);
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