/// @file Решение задачи переноса с неоднородной скоростью.
/// Используется схема высокого порядка (MUSCL) с адаптацией.

#include "fast.h"

// Вектор состояния
struct _U_ {
    double u1, ux, uy, uh, u2;
};

// Для быстрого доступа по типу
_U_ U;

// Векторное поле скорости
Vector3d velocity(const Vector3d& c) {
    return { 1.0, 0.3 + 0.3*std::sin(4 * M_PI * c.x()), 0.0 };
}

// Переменные для сохранения
double get_u(Storage::Item cell)  { return cell(U).u1; }

double get_ux(Storage::Item cell)  { return cell(U).ux; }

double get_uy(Storage::Item cell)  { return cell(U).uy; }

double get_vx(Storage::Item cell) { return velocity(cell.center()).x(); }

double get_vy(Storage::Item cell) { return velocity(cell.center()).y(); }

double compute_dt(ICell& cell) {
    double max_area = 0.0;
    for (auto &face: cell.faces()) {
        max_area = std::max(max_area, face.area());
    }
    double dx = cell.volume() / max_area;
    return dx / velocity(cell.center()).norm();
}

// Посчитаем по Гауссу, но надо бы МНК намутить
// stage in {0, 1}
void compute_grad(ICell& cell, int stage) {
    double ux = 0.0;
    double uy = 0.0;

    double uc = stage < 1 ? cell(U).u1 : cell(U).uh;

    for (auto& face: cell.faces()) {
        auto neib = face.neib();

        double un = stage < 1 ? neib(U).u1 : neib(U).uh;

        Vector3d S = 0.5 * face.normal() * face.area();
        ux += (uc + un) * S.x();
        uy += (uc + un) * S.y();
    }

    cell(U).ux = ux / cell.volume();
    cell(U).uy = uy / cell.volume();
}

double limiter(double t, double b) {
    if (b == 0.0 && t == 0.0) {
        return 0.0;
    }
    double r = t / b;
    return std::max(0.0, std::min(r, 1.0));
}

void fluxes(ICell& cell, double dt, int stage) {
    auto &zc = cell(U);

    double fluxes = 0.0;
    for (auto &face: cell.faces()) {
        auto neib = face.neib();
        auto &zn = neib(U);

        double af = velocity(face.center()).dot(face.normal());
        double a_p = std::max(af, 0.0);
        double a_m = std::min(af, 0.0);

        Vector3d dr = 2.0 * (face.center() - cell.center());

        // Проекции на грань
        double uc = stage < 1 ? zc.u1 : zc.uh;
        double un = stage < 1 ? zn.u1 : zn.uh;

        // Приращения в направлении нормали
        double delta_c = zc.ux * dr.x() + zc.uy * dr.y();
        double delta_n = zn.ux * dr.x() + zn.uy * dr.y();

        double u_m = uc + 0.5 * limiter(delta_c, un - uc) * (un - uc);
        double u_p = un - 0.5 * limiter(un - uc, delta_n) * delta_n;

        fluxes += (a_p * u_m + a_m * u_p) * face.area();
    }

    if (stage < 1) {
        zc.uh = zc.u1 - 0.5 * dt * fluxes / cell.volume();
    }
    else {
        zc.u2 = zc.u1 - dt * fluxes / cell.volume();
    }
}

int main() {
    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables += {"u",  get_u};
    pvd.variables += {"ux", get_ux};
    pvd.variables += {"uy", get_uy};
    pvd.variables += {"vx", get_vx};
    pvd.variables += {"vy", get_vy};

    // Геометрия области
    Rectangle rect(0.0, 1.0, 0.0, 0.6, true);
    rect.set_nx(200);
    rect.set_boundary_flags(
            FaceFlag::PERIODIC, FaceFlag::PERIODIC,
            FaceFlag::PERIODIC, FaceFlag::PERIODIC);

    // Создать сетку
    Mesh mesh(U, &rect);

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
    double CFL = 0.5;

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

        // Определяем dt
        double dt = std::numeric_limits<double>::max();
        for (auto& cell: mesh.cells()) {
            dt = std::min(dt, compute_dt(cell));
        }
        dt *= CFL;

        // Считаем производные
        for (auto& cell: mesh.cells()) {
            compute_grad(cell, 0);
        }

        // Шаг предиктора
        for (auto cell: mesh.cells()) {
            fluxes(cell, dt, 0);
        }

        // Считаем производные
        for (auto& cell: mesh.cells()) {
            compute_grad(cell, 1);
        }

        // Шаг корректора
        for (auto cell: mesh.cells()) {
            fluxes(cell, dt, 1);
        }

        // Обновляем слои
        for (auto cell: mesh.cells()) {
            cell(U).u1 = cell(U).u2;
            cell(U).uh = 0.0;
            cell(U).u2 = 0.0;
        }

        n_step += 1;
        curr_time += dt;
    }

    return 0;
}