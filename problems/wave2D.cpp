/// @file Решение задачи о колебаниях мебраны
/// Используется метод Руенге-Кутта 2 порядка

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
double get_u(Storage::Item cell) {
    return cell(U).u1;
}

double get_bound(Storage::Item cell) {
    if (cell(U).is_bound)
        return 1.0;
    else
        return 0.0;
}

int main() {
    // Файл для записи
    PvdFile pvd("wave2d", R"(C:\Users\Diablo\CLionProjects\zephyr\output\)");

    // Переменные для сохранения
    pvd.variables += {"u", get_u};
    pvd.variables += {"bound", get_bound};

    // Геометрия области
    double h = 0.5;
    Rectangle rect(-1.0, 1.0, -1.0, 1.0, false);
    rect.set_nx(500);
    rect.set_boundary_flags(
            FaceFlag::WALL, FaceFlag::WALL,
            FaceFlag::WALL, FaceFlag::WALL);

    // Создать сетку
    Mesh mesh(U, &rect);

    // Заполняем начальные данные
    for (auto cell: mesh.cells()) {
        auto center = cell.center();
//        cell(U).u1 = (1 - center.x() * center.x() + center.y() * center.y() / r / r) * h;
//        if (abs(center.x()) < 0.2 && abs(center.y()) < 0.2)
//            cell(U).u1 = h * (0.04 - center.x() * center.x()) * (0.04 - center.y() * center.y());
//        else
//            cell(U).u1 = 0.0;
        cell(U).u1 = h * (1 - abs(center.x())) * (1 - abs(center.y()));
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
    double CFL = 0.8;
    double a = 1;

    int n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;
    double finish_time = 1.5;
    double period = 0.05;

    // Определяем dt
    double dt = std::numeric_limits<double>::max();
    for (auto &cell: mesh.cells()) {
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
            pvd.save(mesh.cells(), curr_time);
            next_write += period;
        }

        for (auto cell: mesh.cells()) {
            auto &zc = cell(U);
            if (zc.is_bound) {
                zc.u_av = 0;
                zc.v_av = 0;
                continue;
            }

            double fluxes = 0.0;
            for (auto &face: cell.faces()) {
//                if (face.is_boundary()) {
//                    fluxes -= zc.u1 * face.area();
//                } else {
//                    auto neib = face.neib();
//                    auto &zn = neib(U);
//                    fluxes += (zn.u1 - zc.u1) * face.area() / (cell.center() - face.center()).norm() / 4;
//                }
                auto neib = face.neib();
                auto &zn = neib(U);
                fluxes += (zn.u1 - zc.u1) * face.area() / (cell.center() - face.center()).norm() / 2;
            }

            zc.u_av = zc.u1 + dt / 2 * zc.v1; // u(n + 1/2)
            zc.v_av = zc.v1 + dt / 2 * a * a * fluxes / cell.volume(); // v(n + 1/2)
        }

        for (auto cell: mesh.cells()) {
            auto &zc = cell(U);
            if (zc.is_bound) {
                zc.u2 = 0;
                zc.v2 = 0;
                continue;
            }

            double fluxes = 0.0;
            for (auto &face: cell.faces()) {
//                if (face.is_boundary()) {
//                    fluxes -= zc.u_av * face.area();
//                } else {
//                    auto neib = face.neib();
//                    auto &zn = neib(U);
//                    fluxes += (zn.u_av - zc.u_av) * face.area() / (cell.center() - face.center()).norm() / 4;
//                }
                auto neib = face.neib();
                auto &zn = neib(U);
                fluxes += (zn.u_av - zc.u_av) * face.area() / (cell.center() - face.center()).norm() / 2;
            }

            zc.u2 = zc.u1 + dt * zc.v_av; // u(n + 1)
            zc.v2 = zc.v1 + dt * a * a * fluxes / cell.volume(); // v(n + 1)
        }

        for (auto cell: mesh.cells()) {
            double u = 0, v = 0;
            int count = 0;
            for (auto &face: cell.faces()) {
                count++;
                u += face.neib()(U).u2;
                v += face.neib()(U).v2;
            }
            cell(U).u_av = 0.5 * u / count + 0.5 * cell(U).u2;
            cell(U).v_av = 0.5 * u / count + 0.5 * cell(U).v2;
        }

        // Обновляем слои
        for (auto cell: mesh.cells()) {
            cell(U).u1 = cell(U).u_av;
            cell(U).v1 = cell(U).v_av;
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