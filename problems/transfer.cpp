/// @file transfer.cpp
/// @brief Решение задачи переноса со сложными решателями.

#include <iostream>
#include <iomanip>

#include <zephyr/io/pvd_file.h>
#include <zephyr/geom/vector.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/math/solver/transfer.h>

using namespace zephyr::geom;
using namespace zephyr::mesh;

using zephyr::io::PvdFile;
using zephyr::mesh::generator::Rectangle;

class Solver;

/// @brief Некоторая геметрия
class Body {
public:
    /// @brief Выставить на сетке начальные условия
    virtual void initial(Solver& solver, SoaMesh& mesh, bool exact) const = 0;

    virtual double volume_inside(SoaMesh& body) const = 0;

    std::function<bool(const Vector3d& v)> inside;


    Vector3d C0;  ///< Начальные координаты центра
    Vector3d C;   ///< Текущие координаты центра

    /// @brief Точки границы
    std::vector<Vector3d> vs;
};

/// @brief Начальные условия в виде полосы
class BodyLine : public Body {
public:
    // Полоса задается двумя плоскостями
    Vector3d p1, n1;
    Vector3d p2, n2;


    BodyLine();

    /// @brief Установить начальные данные
    void initial(Solver& solver, SoaMesh& mesh, bool exact) const final;

    double volume_inside(SoaMesh& body) const final;
};

/// @brief Начальные условия в виде квадрата
class BodySquare : public Body {
public:
    double a;     ///< Сторона квадрата

    BodySquare();

    /// @brief Установить начальные данные
    void initial(Solver& solver, SoaMesh& mesh, bool exact) const final;

    double volume_inside(SoaMesh& body) const final;

};

/// @brief Начальные условия в виде круга
class BodyDisk : public Body {
public:
    double R;     ///< Радиус круга

    BodyDisk();

    /// @brief Установить начальные данные
    void initial(Solver& solver, SoaMesh& mesh, bool exact) const final;

    double volume_inside(SoaMesh& body) const final;

};

/// @brief Наследуем собственный решатель от Transfer, теперь переопределив
/// поле скорости можно решать произвольные задачи на перенос.
class Solver : public zephyr::math::Transfer {
public:
    enum class Test {
        Translation,  // Поступательное движение
        Rotation,     // Вращательное движение
    };

    Test test = Test::Translation;

    /// @brief Скорость переноса в соответсвии с Solver::test
    Vector3d velocity(const Vector3d& p) const override;

    /// @brief Сетка с точным решением от времени
    SoaMesh exact(Body& body, double curr_time) const;
};

static Solver::State data;

// Объем тела
double volume(SoaMesh& cells, Storable<double> u1) {
    double sum = 0.0;
    for (auto cell: cells) {
        if (cell.volume() >= 0)
            sum += cell.volume() * cell(u1);
    }
    return sum;
}

int main() {
    // Файл для записи
    PvdFile pvd("mesh", "output");
    PvdFile pvd_body("body", "output");
    PvdFile pvd_exact("exact", "output");

    // Использовать полигональную сетку
    bool voronoi = false;

    // Геометрия области
    Rectangle rect(0.0, 1.0, 0.0, 0.7, voronoi);
    rect.set_nx(200);
    rect.set_boundaries({
        .left   = Boundary::ZOE, .right = Boundary::ZOE,
        .bottom = Boundary::ZOE, .top   = Boundary::ZOE});

    // Создать решатель
    Solver solver;
    solver.set_CFL(0.5);

    // Расщепление по направлениям
    bool splitting = true;

    // Настройки метода
    solver.set_method(Solver::Method::CRP_N1);

    // Настройки теста
    BodyDisk body;
    solver.test = Solver::Test::Translation;

    // Создать сетку
    SoaMesh mesh(rect);

    // Добавить типы
    solver.add_types(mesh);

    data = solver.data;

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables.append("u", data.u1);
    pvd.variables.append("u2", data.u2);
    //pvd.variables.append("n", data.n);
    //pvd.variables.append("p", data.p);
    //pvd.variables += {"du/dx", grad_x};
    //pvd.variables += {"du/dy", grad_y};
    pvd.variables += {"over", [data](QCell& cell) -> double {
        double u = cell(data.u1);
        return u < 0.0 ? u : (u <= 1.0 ? 0.0 / 0.0 : u - 1.0);
    }};
    pvd.variables += {"close", [data](QCell& cell) -> double {
        double u = cell(data.u1);
        return std::abs(u < 0.5 ? u : 1.0 - u);
    }};

    // Настраиваем адаптацию
    mesh.set_max_level(0);
    mesh.set_distributor(solver.distributor());

    // Адаптация под начальные данные
    int n_init_loops = mesh.is_adaptive() ? mesh.max_level() + 2 : 0;
    for (int k = n_init_loops; k >= 0; --k) {
        body.initial(solver, mesh, k < 1);
        solver.set_flags(mesh);
        mesh.refine();
    }

    solver.update_interface(mesh);
    double init_volume = volume(mesh, data.u1);

    int n_step = 0;
    double end_time = 1.0;
    double curr_time = 0.0;
    double write_freq = end_time / 100;
    double write_next = 0.0;

    while (n_step < 1000) {
        if (curr_time >= write_next || curr_time >= end_time) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(8) << std::setprecision(3) << std::fixed
                      << curr_time << ";";

            double curr_volume = volume(mesh, data.u1);
            std::cout << "\tLoss: " << std::setw(12) << std::setprecision(3) << std::scientific
                      << (curr_volume - init_volume) / (init_volume) << "\n";

            pvd.save(mesh, curr_time);

            SoaMesh crop = solver.body(mesh);
            pvd_body.save(crop, curr_time);

            SoaMesh exact = solver.exact(body, curr_time);
            pvd_exact.save(exact, curr_time);

            write_next += write_freq;

            if (curr_time >= end_time) {
                //double vi = body.volume_inside(crop);
                std::cout << std::setprecision(3) << std::fixed;
                //std::cout << "  Volume loss: " << 100 * (1.0 - vi / init_volume) << "%\n";
                //pvd_body.save(crop, curr_time + 1.0e-13);
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
        //solver.set_flags(mesh);

        // Адаптировать сетку
        //mesh.refine();

        n_step += 1;
        curr_time += solver.get_dt();
    }

    return 0;
}

BodyLine::BodyLine() {
    p1 = { 0.1, 0.0, 0.0};
    p2 = { 0.3, 0.0, 0.0};
    n1 = {-1.0, 0.0, 0.0};
    n2 = {+1.0, 0.0, 0.0};

    inside = [this](const Vector3d& v) -> bool {
        return (v - p1).dot(n1) < 0.0 &&
               (v - p2).dot(n2) < 0.0;
    };
}

void BodyLine::initial(Solver& solver, SoaMesh& mesh, bool exact) const {
    for (auto cell: mesh) {
        if (!exact) {
            cell(data.u1) = inside(cell.center());
        }
        else {
            double vol_frac = cell.approx_vol_fraction(inside);
            if (0.0 < vol_frac && vol_frac < 1.0) {
                vol_frac = cell.volume_fraction(inside, 10000);
            }
            cell(data.u1) = vol_frac;
        }
        cell(data.u2) = 0.0;
    }
}

double BodyLine::volume_inside(SoaMesh &body) const {
    double res = 0.0;
    for (auto& cell: body) {
        double a = cell.approx_vol_fraction(inside);
        if (0.0 < a && a < 1.0) {
            a = cell.volume_fraction(inside, 1000);
        }
        res += a * cell.volume();
    }
    return res;
}

BodySquare::BodySquare() {
    a = 0.2;
    C = C0 = {0.15, 0.5, 0.0};

    // характеристическая функция области
    inside = [this](const Vector3d &v) -> bool {
        return std::abs(v.x() - C.x()) < 0.5 * a &&
               std::abs(v.y() - C.y()) < 0.5 * a;
    };

    // граница области
    vs = {
            Vector3d{C.x() - 0.5 * a, C.y() - 0.5 * a, 0.0},
            Vector3d{C.x() + 0.5 * a, C.y() - 0.5 * a, 0.0},
            Vector3d{C.x() + 0.5 * a, C.y() + 0.5 * a, 0.0},
            Vector3d{C.x() - 0.5 * a, C.y() + 0.5 * a, 0.0},
    };
}

void BodySquare::initial(Solver& solver, SoaMesh& mesh, bool exact) const {
    for (auto cell: mesh) {
        if (!exact) {
            cell(data.u1) = inside(cell.center());
        }
        else {
            double vol_frac = cell.approx_vol_fraction(inside);
            if (0.0 < vol_frac && vol_frac < 1.0) {
                vol_frac = cell.volume_fraction(inside, 10000);
            }
            cell(data.u1) = vol_frac;
        }
        cell(data.u2) = 0.0;
    }
}

double BodySquare::volume_inside(SoaMesh &body) const {
    double res = 0.0;
    for (auto& cell: body) {
        double vol = cell.approx_vol_fraction(inside);
        if (0.0 < vol && vol < 1.0) {
            vol = cell.volume_fraction(inside, 1000);
        }
        cell(data.u1) = vol;
        vol *= cell.volume();
        res += vol;
    }
    return res;
}

BodyDisk::BodyDisk() {
    R = 0.1;
    C = C0 = {0.15, 0.5, 0.0};

    size_t M = 100;
    vs.resize(M);
    for (size_t i = 0; i < M; ++i) {
        double phi = 2.0 * M_PI * (i - 1.0) / M;
        vs[i].x() = C.x() + R * std::cos(phi);
        vs[i].y() = C.y() + R * std::sin(phi);
        vs[i].z() = 0.0;
    }
}

void BodyDisk::initial(Solver& solver, SoaMesh& mesh, bool exact) const {
    for (auto cell: mesh) {
        auto poly = cell.polygon();
        cell(data.u1) = poly.disk_clip_area(C, R) / cell.volume();
        cell(data.u2) = 0.0;
    }
}

double BodyDisk::volume_inside(SoaMesh &body) const {
    double res = 0.0;
    for (auto& cell: body) {
        auto poly = cell.polygon();
        double vol = poly.disk_clip_area(C, R);
        cell(data.u1) = vol / cell.volume();
        res += vol;
    }
    return res;
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

SoaMesh Solver::exact(Body& body, double curr_time) const {
    using zephyr::geom::Quad;
    using zephyr::mesh::AmrCell;

    // Точная граница в начальный момент времени
    std::vector<Vector3d> vs = body.vs;

    // Точное решение
    if (test == Test::Translation) {
        Vector3d V0 = {0.7, -0.35, 0.0};
        for (auto& v: vs) {
            v += curr_time * V0;
        }
        body.C = body.C0 + curr_time * V0;
    }
    else if (test == Test::Rotation) {
        Vector3d center = {0.5, 0.5, 0.0};  // Центр вращения
        Vector3d omega = {0.0, 0.0, M_PI};  // Угловая частота

        double phi = -M_PI * curr_time;
        for (auto& v: vs) {
            Vector3d r = v - center;
            v.x() = center.x() + r.x() * std::cos(phi) + r.y() * std::sin(phi);
            v.y() = center.y() - r.x() * std::sin(phi) + r.y() * std::cos(phi);
            v.z() = 0.0;
        }

        Vector3d r = body.C0 - center;
        body.C.x() = center.x() + r.x() * std::cos(phi) + r.y() * std::sin(phi);
        body.C.y() = center.y() - r.x() * std::sin(phi) + r.y() * std::cos(phi);
        body.C.z() = 0.0;
    }

    SoaMesh cells(2, false);
    for (size_t i = 0; i < vs.size(); ++i) {
        size_t j = (i + 1) % vs.size();
        Line line = {vs[i], vs[j]};
        cells.push_back(line);
    }
    return cells;
}