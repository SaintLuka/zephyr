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
        return { 1.0, 0.3 + 0.3*std::sin(4 * M_PI * c.x()), 0.0 };
    }
};

// Получим тип данных
Solver::State U = Solver::datatype();

// Переменные для сохранения
double get_u(Storage::Item cell)  { return cell(U).u1; }

double get_ux(Storage::Item cell)  { return cell(U).ux; }

double get_uy(Storage::Item cell)  { return cell(U).uy; }

int main() {
    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables += {"u",  get_u};
    pvd.variables += {"ux", get_ux};
    pvd.variables += {"uy", get_uy};

    // Геометрия области
    Rectangle rect(0.0, 1.0, 0.0, 0.6, true);
    rect.set_nx(200);
    rect.set_boundary_flags(
            FaceFlag::PERIODIC, FaceFlag::PERIODIC,
            FaceFlag::PERIODIC, FaceFlag::PERIODIC);

    // Создать сетку
    Mesh mesh(U, &rect);

    // Создать решатель
    Solver solver;

    // Настроим решатель
    solver.set_CFL(0.5);
    solver.set_accuracy(3);
    solver.set_limiter("van Leer");

    // Заполняем начальные данные
    Vector3d v1(rect.x_min(), rect.y_min(), 0.0);
    Vector3d v2(rect.x_max(), rect.y_max(), 0.0);
    Vector3d vc = 0.5 * (v1 + v2);
    double D = 0.1 * (v2 - v1).norm();
    for (auto cell: mesh.cells()) {
        cell(U).u1 = (cell.center() - vc).norm() < D ? 1.0 : 0.0;
        cell(U).u2 = 0.0;
    }

    // Число Куранта
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

        solver.update(mesh);

        n_step += 1;
        curr_time += solver.dt();
    }

    return 0;
}