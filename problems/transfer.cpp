/// @file Решение задачи переноса с CRP решателем.

#include "fast.h"

#include <zephyr/math/solver/transfer.h>

using namespace zephyr::geom;

/// @brief Наследуем собственный решатель от Transfer, теперь переопределив
/// поле скорости можно решать произвольные задачи на перенос.
class Solver : public zephyr::math::Transfer {
public:

    // Параметры для тестирования
    enum class Init {
        Line,
        Square,
        Disk,
    };

    enum class Test {
        Translation,  // Поступательное движение
        Rotation,     // Вращательное движение
    };

    Init init = Init::Square;
    Test test = Test::Translation;

    /// @brief Установить начальные данные в соответствии
    /// с указаной в Solver::init
    void setup_initial(EuMesh& mesh, bool exact = false);

    /// @brief Скорость переноса в соответсвии с Solver::test
    Vector3d velocity(const Vector3d& p) const override;

    /// @brief Сетка с точным решением от времени
    AmrStorage exact(double curr_time) const;

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
    PvdFile pvd_exact("exact", "output");
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

    // Использовать полигональную сетку
    bool voronoi = false;

    // Геометрия области
    Rectangle rect(0.0, 1.0, 0.0, 0.7, voronoi);
    rect.set_nx(100);
    rect.set_boundaries({
        .left   = Boundary::ZOE, .right = Boundary::ZOE,
        .bottom = Boundary::ZOE, .top   = Boundary::ZOE});

    // Создать решатель
    Solver solver;
    solver.set_CFL(0.5);

    // Расщепление по направлениям
    bool splitting = false;

    // Настройки метода
    solver.set_method(Solver::Method::CRP_N);

    // Настройки теста
    solver.init = Solver::Init::Disk;
    solver.test = Solver::Test::Rotation;

    // Создать сетку
    EuMesh mesh(U, &rect);

    // Настраиваем адаптацию
    mesh.set_max_level(0);
    mesh.set_distributor(solver.distributor());

    // Адаптация под начальные данные
    int n_init_loops = mesh.is_adaptive() ? mesh.max_level() + 2 : 0;
    for (int k = n_init_loops; k >= 0; --k) {
        solver.setup_initial(mesh, k < 1);
        solver.set_flags(mesh);
        mesh.refine();
    }
    solver.update_interface(mesh);

    double init_volume = volume(mesh);
    std::cout << "Начальный объем: " << init_volume << "\n";

    int n_step = 0;
    double end_time = 1.0;
    double curr_time = 0.0;
    double write_freq = end_time / 100;
    double write_next = 0.0;

    while (curr_time <= end_time && n_step < 1000) {
        if (curr_time >= write_next || curr_time >= end_time) {
            std::cout << "\tШаг: " << std::setw(6) << n_step << ";"
                      << "\tВремя: " << std::setw(8) << std::setprecision(3) << std::fixed
                      << curr_time << ";";

            double curr_volume = volume(mesh);
            std::cout << "\tОшибка: " << std::setw(12) << std::setprecision(3) << std::scientific
                      << (curr_volume - init_volume) / (init_volume) << "\n";

            pvd.save(mesh.locals(), curr_time);

            AmrStorage body = solver.body(mesh);
            pvd_body.save(body, curr_time);

            AmrStorage exact = solver.exact(curr_time);
            pvd_exact.save(exact, curr_time);

            //AmrStorage scheme = solver.scheme(mesh);
            //pvd_scheme.save(scheme, curr_time);

            write_next += write_freq;

            if (curr_time >= end_time) {
                break;
            }
        }

        // Определить шаг
        double dt = solver.compute_dt(mesh);
        if (curr_time + dt > end_time) {
            dt = end_time - curr_time;
        }
        solver.set_dt(dt);

        // Шаг интегрирования
        if (splitting) {
            solver.update(mesh, Direction::X);
            solver.update(mesh, Direction::Y);
        }
        else {
            solver.update(mesh, Direction::ANY);
        }

        // Установить флаги адаптации
        solver.set_flags(mesh);

        // Адаптировать сетку
        mesh.refine();

        n_step += 1;
        curr_time += solver.get_dt();
    }

    return 0;
}

// Начальное условие в виде полосы
void setup_initial_line(EuMesh& mesh, bool exact = false) {
    // Полоса задана двумя линиями
    auto inside = [](const Vector3d& v) -> bool {
        const Vector3d p1 = {0.1, 0.0, 0.0};
        const Vector3d p2 = {0.3, 0.0, 0.0};
        const Vector3d n1 = {-1.0, 0.0, 0.0};
        const Vector3d n2 = {+1.0, 0.0, 0.0};

        return (v - p1).dot(n1) < 0.0 &&
               (v - p2).dot(n2) < 0.0;
    };

    for (auto cell: mesh) {
        if (!exact) {
            cell(U).u1 = inside(cell.center());
        }
        else {
            double vol_frac = cell.approx_vol_fraction(inside);
            if (0.0 < vol_frac && vol_frac < 1.0) {
                vol_frac = cell.volume_fraction(inside, 10000);
            }
            cell(U).u1 = vol_frac;
        }
        cell(U).u2 = 0.0;
    }
}

// Начальное условие в виде квадрата
void setup_initial_square(EuMesh& mesh, bool exact = false) {
    auto inside = [](const Vector3d& v) -> bool {
        const double a = 0.2;
        const double xc = 0.15;
        const double yc = 0.5;

        return std::abs(v.x() - xc) < 0.5 * a &&
               std::abs(v.y() - yc) < 0.5 * a;
    };

    for (auto cell: mesh) {
        if (!exact) {
            cell(U).u1 = inside(cell.center());
        }
        else {
            double vol_frac = cell.approx_vol_fraction(inside);
            if (0.0 < vol_frac && vol_frac < 1.0) {
                vol_frac = cell.volume_fraction(inside, 10000);
            }
            cell(U).u1 = vol_frac;
        }
        cell(U).u2 = 0.0;
    }
}

// Начальное условие в виде круга
void setup_initial_disk(EuMesh& mesh, bool exact = true) {
    const double R = 0.1;
    const Vector3d c = {0.15, 0.5, 0.0};

    for (auto cell: mesh) {
        auto poly = cell.polygon();
        cell(U).u1 = poly.disk_clip_area(c, R) / cell.volume();
        cell(U).u2 = 0.0;
    }
}

void Solver::setup_initial(EuMesh& mesh, bool exact) {
    switch (init) {
        case Init::Line:
            setup_initial_line(mesh, exact);
            return;
        case Init::Square:
            setup_initial_square(mesh, exact);
            return;
        case Init::Disk:
            setup_initial_disk(mesh, exact);
            return;
    }
}

Vector3d Solver::velocity(const Vector3d& p) const {
    if (test == Test::Translation) {
        Vector3d V0 = {0.7, -0.35, 0.0};
        return V0;
    }
    else if (test == Test::Rotation) {
        Vector3d center = {0.5, 0.5, 0.0};  // Центр вращения
        Vector3d omega = {0.0, 0.0, M_PI};  // Угловая частота
        return omega.cross(p - center);
    }
    else {
        return Vector3d::Zero();
    }
}

AmrStorage Solver::exact(double curr_time) const {
    using zephyr::geom::Quad;
    using zephyr::geom::AmrCell;

    // Точная граница
    std::vector<Vector3d> vs;

    // Установим начальные данные
    if (init == Init::Disk) {
        double R = 0.1;
        Vector3d c = {0.15, 0.5, 0.0};

        size_t M = 100;
        vs.resize(M);
        for (size_t i = 0; i < M; ++i) {
            double phi = 2.0 * M_PI * (i - 1.0) / M;
            vs[i].x() = c.x() + R * std::cos(phi);
            vs[i].y() = c.y() + R * std::sin(phi);
            vs[i].z() = 0.0;
        }
    }
    else if (init == Init::Square) {
        const double a = 0.2;
        const double xc = 0.15;
        const double yc = 0.5;

        vs = {
                Vector3d{xc - 0.5 * a, yc - 0.5 * a, 0.0},
                Vector3d{xc + 0.5 * a, yc - 0.5 * a, 0.0},
                Vector3d{xc + 0.5 * a, yc + 0.5 * a, 0.0},
                Vector3d{xc - 0.5 * a, yc + 0.5 * a, 0.0},
        };
    }

    // Точное решение
    if (test == Test::Translation) {
        Vector3d V0 = {0.7, -0.35, 0.0};
        for (auto& v: vs) {
            v += curr_time * V0;
        }
    }
    else if (test == Test::Rotation) {
        Vector3d center = {0.5, 0.5, 0.0};  // Центр вращения
        Vector3d omega = {0.0, 0.0, M_PI};  // Угловая частота
        for (auto& v: vs) {
            Vector3d r = v - center;
            double phi = -M_PI * curr_time;
            v.x() = center.x() + r.x() * std::cos(phi) + r.y() * std::sin(phi);
            v.y() = center.y() - r.x() * std::sin(phi) + r.y() * std::cos(phi);
            v.z() = 0.0;
        }
    }

    AmrStorage cells(vs.size());
    for (size_t i = 0; i < vs.size(); ++i) {
        size_t j = (i + 1) % vs.size();
        Quad quad = {
                vs[i], vs[j], vs[i], vs[j]
        };
        AmrCell cell(quad);
        cells[i] = cell;
    }

    return cells;

}