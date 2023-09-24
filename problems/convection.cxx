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
        return { 0.8, 0.40 + 0.3*std::sin(4 * M_PI * c.x()), 0.0 };
    }
};

// Получим тип данных
Solver::State U = Solver::datatype();

// Переменные для сохранения
double get_u(Storage::Item cell)  { return cell(U).u1; }

double get_ux(Storage::Item cell)  { return cell(U).ux; }

double get_uy(Storage::Item cell)  { return cell(U).uy; }

double get_lvl(const Storage::Item &cell)  { return cell.level(); }

const double margin = 0.0198;

// Начальное условие в виде круга
void setup_initial_1(Mesh& mesh, double D) {
    double R = D / 2.0;
    Vector3d vc = {R + margin, R + margin, 0.0};
    for (auto cell: mesh.cells()) {
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

    for (auto cell: mesh.cells()) {
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
    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables += {"u",  get_u};
    pvd.variables += {"ux", get_ux};
    pvd.variables += {"uy", get_uy};
    pvd.variables += {"lvl", get_lvl};

    // Геометрия области
    Rectangle rect(0.0, 1.0, 0.0, 0.6, false);
    rect.set_nx(50);
    rect.set_boundary_flags(
            FaceFlag::ZOE, FaceFlag::ZOE,
            FaceFlag::ZOE, FaceFlag::ZOE);

    // Создать решатель
    Solver solver;

    // Настроим решатель
    solver.set_CFL(0.5);
    solver.set_accuracy(3);
    solver.set_limiter("van Leer");

    // Создать сетку
    Mesh mesh(U, &rect);

    // Настраиваем адаптацию
    mesh.set_max_level(5);
    mesh.set_distributor(solver.distributor());

    // Заполняем начальные данные
    Vector3d v_min(rect.x_min(), rect.y_min(), 0.0);
    Vector3d v_max(rect.x_max(), rect.y_max(), 0.0);
    double D = 0.1*(v_max - v_min).norm();

    // Адаптация под начальные данные
    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        setup_initial_1(mesh, D);
        solver.set_flags(mesh);
        mesh.refine();
    }

    int n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    while(curr_time < 1.0) {
        if (curr_time >= next_write) {
            std::cout << "\tШаг: " << std::setw(6) << n_step << ";"
                      << "\tВремя: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh.cells(), curr_time);
            next_write += 0.02;
        }

        // Шаг решения
        solver.update(mesh);

        // Установить флаги адаптации
        solver.set_flags(mesh);

        // Адаптировать сетку
        mesh.refine();

        n_step += 1;
        curr_time += solver.dt();
    }

    return 0;
}