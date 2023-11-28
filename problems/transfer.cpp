/// @file Решение задачи переноса с CRP решателем.

#include "fast.h"

#include <zephyr/math/solver/transfer.h>

/// @brief Наследуем собственный решатель от Transfer, теперь переопределив
/// поле скорости можно решать произвольные задачи на перенос.
class Solver : public zephyr::math::Transfer {
public:

    /// @brief Задаем скорость переноса
    Vector3d velocity(const Vector3d& c) const override {
        return {0.5, 0.5, 0.0};

        double rx = c.x() - 0.5;
        double ry = c.y() - 0.5;
        double phi = std::atan2(ry, rx);
        double v = 0.5 * std::sqrt(rx * rx + ry * ry);
        return {v * std::sin(phi), -v * std::cos(phi), 0.0};

         //0.40 + 0.3*std::sin(4 * M_PI * c.x()), 0.0 };
    }
};

// Получим тип данных
Solver::State U = Solver::datatype();

// Переменные для сохранения
double get_u(AmrStorage::Item& cell)  { return cell(U).u1; }

double get_lvl(AmrStorage::Item& cell)  { return cell.level; }

double normal_x(AmrStorage::Item& cell) { return cell(U).n.x(); }

double normal_y(AmrStorage::Item& cell) { return cell(U).n.y(); }

double point_x(AmrStorage::Item& cell) { return cell(U).p.x(); }

double point_y(AmrStorage::Item& cell) { return cell(U).p.y(); }

double get_over(AmrStorage::Item& cell) {
    auto u = cell(U).u1;
    if (u < 0.0) {
        return u;
    }
    else if (u <= 1.0) {
        return 0.0 / 0.0;
    }
    else {
        return u - 1.0;
    }
}

double get_close(AmrStorage::Item& cell) {
    auto u = cell(U).u1;
    return std::abs(u < 0.5 ? u : 1.0 - u);
}

const double margin = 0.2198;

inline double sqr(double x) {
    return x * x;
}

// Начальное условие в виде полосы
void setup_initial_0(EuMesh& mesh, double D) {
    // Обычной блок с резкими границами
    auto func1 = [](double x) -> double {
        return (0.1 < x && x < 0.2) ? 1.0 : 0.0;
    };
    // Блок с гладкими границами
    auto func2 = [](double x) -> double {
        if (x < 0.1) {
            return 0.0;
        }
        else if (x < 0.2) {
            return sqr(std::sin(0.5 * M_PI * (x - 0.1) / 0.1));
        }
        else if (x < 0.4) {
            return 1.0;
        }
        else if (x < 0.5) {
            return sqr(std::sin(0.5 * M_PI * (0.5 - x) / 0.1));
        }
        else {
            return 0.0;
        }
    };
    // Блок с гладкими границами
    auto func3 = [](double x) -> double {
        if (x < 0.1) {
            return 0.0;
        }
        else if (x < 0.5) {
            return 0.4 + 0.6*std::pow(1.0 - std::cos(1.0 * M_PI * (0.3 - x) / 0.4), 0.4);
        }
        else {
            return 0.0;
        }
    };

    for (auto cell: mesh) {
        double x = cell.center().x();
        cell(U).u1 = func2(x);
        cell(U).u2 = 0.0;
    }
}

// Начальное условие в виде круга
void setup_initial_1(EuMesh& mesh, double D) {
    double R = D / 2.0;
    Vector3d vc = {R + margin, R + margin, 0.0};
    for (auto cell: mesh) {
        cell(U).u1 = (cell.center() - vc).norm() < R ? 1.0 : 0.0;
        cell(U).u2 = 0.0;
    }
}

// Начальное условие в виде квадрата
void setup_initial_2(EuMesh& mesh, double D) {
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

// Объем тела
double volume(EuMesh& cells) {
    double sum = 0.0;
    for (auto cell: cells) {
        sum += cell.volume() * cell(U).u1;
    }
    return sum;
}

int main() {
    // Файл для записи
    PvdFile pvd("mesh", "output");
    PvdFile pvd_body("body", "output");
    PvdFile pvd_scheme("scheme", "output");

    // Переменные для сохранения
    pvd.variables += {"u",  get_u};
    pvd.variables += {"lvl", get_lvl};
    pvd.variables += {"n.x", normal_x};
    pvd.variables += {"n.y", normal_y};
    pvd.variables += {"p.x", point_x};
    pvd.variables += {"p.y", point_y};
    pvd.variables += {"over", get_over};
    pvd.variables += {"close", get_close};

    pvd_scheme.variables += {"u",  get_u};

    // Геометрия области
    Rectangle rect(0.0, 1.0, 0.0, 1.0, false);
    rect.set_nx(200);
    rect.set_boundaries({
        .left   = Boundary::ZOE, .right = Boundary::ZOE,
        .bottom = Boundary::ZOE, .top   = Boundary::ZOE});

    // Создать решатель
    Solver solver;
    solver.set_CFL(0.5);
    solver.set_version(1);
    solver.dir_splitting(true);

    // Создать сетку
    EuMesh mesh(U, &rect);

    // Настраиваем адаптацию
    mesh.set_max_level(0);
    mesh.set_distributor(solver.distributor());

    // Заполняем начальные данные
    Vector3d v_min(rect.x_min(), rect.y_min(), 0.0);
    Vector3d v_max(rect.x_max(), rect.y_max(), 0.0);
    double D = 0.2*(v_max - v_min).norm();

    // Адаптация под начальные данные
    for (int k = 0; k < mesh.max_level() + 3; ++k) {
        setup_initial_2(mesh, D);
        solver.set_flags(mesh);
        mesh.refine();
    }

    solver.compute_normals(mesh);
    solver.find_sections(mesh);

    double init_volume = volume(mesh);
    std::cout << "Начальный объем: " << init_volume << "\n";

    int n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    while(curr_time < 10.0 && n_step < 200) {
        if (curr_time >= next_write) {
            solver.compute_normals(mesh);
            solver.find_sections(mesh);

            std::cout << "\tШаг: " << std::setw(6) << n_step << ";"
                      << "\tВремя: " << std::setw(8) << std::setprecision(3) << std::fixed
                      << curr_time << ";";

            double curr_volume = volume(mesh);
            std::cout << "\tОшибка: " << std::setw(12) << std::setprecision(3) << std::scientific
                      << (curr_volume - init_volume) / (init_volume) << "\n";

            pvd.save(mesh.locals(), curr_time);

            AmrStorage body = solver.body(mesh);
            pvd_body.save(body, curr_time);

            AmrStorage scheme = solver.scheme(mesh);
            pvd_scheme.save(scheme, curr_time);

            next_write += 0.0;
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