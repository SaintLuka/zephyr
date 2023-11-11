/// @file Решение задачи переноса с неоднородной скоростью.
/// Используется простейшая схема первого порядка (upwind)

#include "fast.h"

// Вектор состояния
struct _U_ {
    double u1, u2;
};

// Для быстрого доступа по типу
_U_ U;

// Векторное поле скорости
Vector3d velocity(const Vector3d& c) {
    return { 1.0, 0.3 + 0.3*std::sin(4 * M_PI * c.x()), 0.0 };
}

// Переменные для сохранения
double get_u(AmrStorage::Item& cell)  {
    return cell(U).u1;
}

double get_vx(AmrStorage::Item& cell) {
    return velocity(cell.center).x();
}

double get_vy(AmrStorage::Item& cell) {
    return velocity(cell.center).y();
}

int main() {
    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables += {"u",  get_u};
    pvd.variables += {"vx", get_vx};
    pvd.variables += {"vy", get_vy};

    // Геометрия области
    Rectangle rect(0.0, 1.0, 0.0, 0.6, true);
    rect.set_nx(200);
    rect.set_boundaries({
            .left   = Boundary::PERIODIC, .right = Boundary::PERIODIC,
            .bottom = Boundary::PERIODIC, .top   = Boundary::PERIODIC});

    // Создать сетку
    EuMesh mesh(U, &rect);

    // Заполняем начальные данные
    Vector3d v1(rect.x_min(), rect.y_min(), 0.0);
    Vector3d v2(rect.x_max(), rect.y_max(), 0.0);
    Vector3d vc = 0.5 * (v1 + v2);
    double D = 0.1 * (v2 - v1).norm();
    for (auto cell: mesh) {
        cell(U).u1 = (cell.center() - vc).norm() < D ? 1.0 : 0.0;
        cell(U).u2 = 0.0;
    }

    // Число Куранта
    double CFL = 0.5;

    int n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    while(curr_time <= 1.0) {
        if (curr_time >= next_write) {
            std::cout << "\tШаг: " << std::setw(6) << n_step << ";"
                      << "\tВремя: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh.locals(), curr_time);
            next_write += 0.02;
        }

        // Определяем dt
        double dt = std::numeric_limits<double>::max();
        for (auto& cell: mesh) {
            double max_area = 0.0;
            for (auto &face: cell.faces()) {
                max_area = std::max(max_area, face.area());
            }
            double dx = cell.volume() / max_area;
            dt = std::min(dt, dx / velocity(cell.center()).norm());
        }
        dt *= CFL;

        // Расчет по схеме upwind
        for (auto cell: mesh) {
            auto& zc = cell(U);

            double fluxes = 0.0;
            for (auto& face: cell.faces()) {
                auto neib = face.neib();
                auto& zn = neib(U);

                double af = velocity(face.center()).dot(face.normal());
                double a_p = std::max(af, 0.0);
                double a_m = std::min(af, 0.0);

                fluxes += (a_p * zc.u1 + a_m * zn.u1) * face.area();
            }

            zc.u2 = zc.u1 - dt * fluxes / cell.volume();
        }

        // Обновляем слои
        for (auto cell: mesh) {
            cell(U).u1 = cell(U).u2;
            cell(U).u2 = 0.0;
        }

        n_step += 1;
        curr_time += dt;
    }

    return 0;
}