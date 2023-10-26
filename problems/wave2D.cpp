/// @file Решение задачи о колебаниях мебраны
/// Используется метод Рунге-Кутта 2 порядка

#include "fast.h"


// Вектор состояния
struct _U_ {
    double u1 = 0.0, u_av = 0.0, u2 = 0.0;
    double v1 = 0.0, v_av = 0.0, v2 = 0.0;
    bool is_bound = false;
};

// Для быстрого доступа по типу
_U_ U;

// Переменные для сохранения
double get_u(AmrStorage::Item& cell) {
    return cell(U).u1;
}

double get_bound(AmrStorage::Item& cell) {
    if (cell(U).is_bound)
        return 1.0;
    else
        return 0.0;
}

int main() {
    // Файл для записи
    PvdFile pvd("wave2d", "output");

    // Переменные для сохранения
    pvd.variables += {"u", get_u};
    pvd.variables += {"bound", get_bound};

    // Геометрия области
    Rectangle rect(-1.0, 1.0, -1.0, 1.0, false);
    rect.set_nx(200);
    rect.set_boundaries({
        .left   = Boundary::WALL, .right = Boundary::WALL,
        .bottom = Boundary::WALL, .top   = Boundary::WALL});

    // Создать сетку
    Mesh mesh(U, &rect);

    // Заполняем начальные данные
    Vector3d v1(rect.x_min(), rect.y_min(), 0.0);
    Vector3d v2(rect.x_max(), rect.y_max(), 0.0);
    Vector3d vc = 0.5 * (v1 + v2);
    double D = 0.1 * (v2 - v1).norm();
    for (auto cell: mesh) {
        auto center = cell.center();
        double r = (center - vc).norm() / D;
        double cos = std::cos(M_PI * r / 2);
        cell(U).u1 = r < 1.0 ? cos*cos*cos : 0.0;
        cell(U).u2 = 0.0;
        cell(U).v1 = 0.0;
        cell(U).v2 = 0.0;
        for (auto &face: cell.faces()) {
            if (face.is_boundary()) {
                cell(U).is_bound = true;
                cell(U).u1 = 0;
                break;
            }
        }
    }

    // Число Куранта
    double CFL = 0.3;
    double a = 1;

    int n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;
    double finish_time = 1.5;
    double period = 0.05;

    // Определяем dt
    double dt = std::numeric_limits<double>::max();
    for (auto &cell: mesh) {
        double max_area = 0.0;
        for (auto &face: cell.faces()) {
            max_area = std::max(max_area, face.area());
        }
        dt = std::min(dt, cell.volume() / max_area);
    }
    dt *= CFL / a / 2;

    while (curr_time <= finish_time) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh.locals(), curr_time);
            next_write += period;
        }

        for (auto cell: mesh) {
            auto &zc = cell(U);
            if (zc.is_bound) { 
                zc.u_av = 0;
                zc.v_av = 0;
                continue;
            }

            double fluxes = 0.0;
            for (auto &face: cell.faces()) {
                auto neib = face.neib();
                auto &zn = neib(U);
                fluxes += 0.5 * (zn.u1 - zc.u1) * face.area() / (cell.center() - face.center()).norm();
            }

            zc.u_av = zc.u1 + 0.5 * dt * zc.v1; // u(n + 1/2)
            zc.v_av = zc.v1 + 0.5 * dt * a * a * fluxes / cell.volume(); // v(n + 1/2)
        }

        for (auto cell: mesh) {
            auto &zc = cell(U);
            if (zc.is_bound) {
                zc.u2 = 0;
                zc.v2 = 0;
                continue;
            }

            double fluxes = 0.0;
            for (auto &face: cell.faces()) {
                auto neib = face.neib();
                auto &zn = neib(U);
                fluxes += 0.5 * (zn.u_av - zc.u_av) * face.area() / (cell.center() - face.center()).norm();
            }

            zc.u2 = zc.u1 + dt * zc.v_av; // u(n + 1)
            zc.v2 = zc.v1 + dt * a * a * fluxes / cell.volume(); // v(n + 1)
        }

        // Обновляем слои
        for (auto cell: mesh) {
            cell(U).u1 = cell(U).u2;
            cell(U).v1 = cell(U).v2;
            cell(U).u_av = 0.0;
            cell(U).v_av = 0.0;
            cell(U).u2 = 0.0;
            cell(U).v2 = 0.0;
        }

        n_step += 1;
        curr_time += dt;
    }

    return 0;
}