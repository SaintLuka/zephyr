/// @file Решение задачи переноса с CRP решателем.

#include <iostream>
#include <iomanip>

#include <zephyr/geom/vector.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/io/pvd_file.h>
#include <zephyr/math/solver/transfer.h>

using namespace zephyr::geom;

using zephyr::io::PvdFile;
using zephyr::mesh::generator::Rectangle;

class Solver;

/// @brief Некоторая геметрия
class Body {
public:
    /// @brief Выставить на сетке начальные условия
    virtual void initial(Solver& solver, EuMesh& mesh, bool exact) const = 0;

    virtual double volume_inside(AmrStorage& body) const = 0;

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
    void initial(Solver& solver, EuMesh& mesh, bool exact) const final;

    double volume_inside(AmrStorage& body) const final;
};

/// @brief Начальные условия в виде квадрата
class BodySquare : public Body {
public:
    double a;     ///< Сторона квадрата

    BodySquare();

    /// @brief Установить начальные данные
    void initial(Solver& solver, EuMesh& mesh, bool exact) const final;

    double volume_inside(AmrStorage& body) const final;

};

/// @brief Начальные условия в виде круга
class BodyDisk : public Body {
public:
    double R;     ///< Радиус круга

    BodyDisk();

    /// @brief Установить начальные данные
    void initial(Solver& solver, EuMesh& mesh, bool exact) const final;

    double volume_inside(AmrStorage& body) const final;

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
    AmrStorage exact(Body& body, double curr_time) const;
};

// Получим тип данных
Solver::State U = Solver::datatype();

// Переменные для сохранения
double get_u(AmrStorage::Item& cell)  { return cell(U).u1; }

double get_lvl(AmrStorage::Item& cell)  { return cell.level; }

double grad_x(AmrStorage::Item& cell) { return cell(U).du_dx; }

double grad_y(AmrStorage::Item& cell) { return cell(U).du_dy; }

double normal_x(AmrStorage::Item& cell) { return cell(U).n.x(); }

double normal_y(AmrStorage::Item& cell) { return cell(U).n.y(); }

double point_x(AmrStorage::Item& cell) { return cell(U).p.x(); }

double point_y(AmrStorage::Item& cell) { return cell(U).p.y(); }

double get_over(AmrStorage::Item& cell) {
    auto u = cell(U).u1;
    return u < 0.0 ? u : (u <= 1.0 ? 0.0 / 0.0 : u - 1.0);
}
double get_fr(AmrStorage::Item& cell){
    return  cell(U).flow_u[0];
}
double get_ft(AmrStorage::Item& cell){
    return  cell(U).flow_u[1];
}
double get_fl(AmrStorage::Item& cell){
    return  cell(U).flow_u[2];
}
double get_fb(AmrStorage::Item& cell){
    return  cell(U).flow_u[3];
}
double get_u2(AmrStorage::Item& cell){
    return  cell(U).u2;
}
double get_u3(AmrStorage::Item& cell){
    return  (cell(U).flow_u[3] + cell(U).flow_u[2]  + cell(U).flow_u[1] + cell(U).flow_u[0] +  cell(U).u1)/5 ;
}
double get_close(AmrStorage::Item& cell) {
    auto u = cell(U).u1;
    return std::abs(u < 0.5 ? u : 1.0 - u);
}
double get_lambda_alpha(AmrStorage::Item& cell) {
    return cell(U).lambda_alpha;
}
double get_lambda_beta(AmrStorage::Item& cell) {
    return cell(U).lambda_beta;
}
double get_sp(AmrStorage::Item& cell){
    if (cell(U).lambda_beta*cell(U).lambda_alpha >0) return 1; // 1 - usual point
    return 0; // sonic point
}

// Объем тела
double volume(EuMesh& cells) {
    double sum = 0.0;
    for (auto cell: cells) {
        if (cell.volume() >= 0)
            sum += cell.volume() * cell(U).u1;
    }
    return sum;
}

int main() {
    // Файл для записи
    PvdFile pvd("mesh", "output");
    PvdFile pvd_body("body", "output");
    PvdFile pvd_exact("exact", "output");
    //PvdFile pvd_scheme("scheme", "output");

    // Переменные для сохранения
    pvd.variables += {"u",  get_u};
    pvd.variables += {"u2",  get_u2};
    pvd.variables += {"u3",  get_u3};
    /*
    pvd.variables += {"fr",  get_fr};
    pvd.variables += {"ft",  get_ft};
    pvd.variables += {"fl",  get_fl};
    pvd.variables += {"fb",  get_fb};
     */
    pvd.variables += {"lvl", get_lvl};
    pvd.variables += {"n.x", normal_x};
    pvd.variables += {"n.y", normal_y};
    //pvd.variables += {"du/dx", grad_x};
    //pvd.variables += {"du/dy", grad_y};
    pvd.variables += {"p.x", point_x};
    pvd.variables += {"p.y", point_y};
    pvd.variables += {"over", get_over};
    pvd.variables += {"close", get_close};

    pvd_body.variables = pvd.variables;

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
    solver.set_CFL(0.13);
    solver.set_mnt(true);

    // Расщепление по направлениям
    bool splitting = true;

    // Настройки метода
    solver.set_method(Solver::Method::CRP_N1);

    // Настройки теста
    BodySquare body;
    solver.test = Solver::Test::Translation;

    // Создать сетку
    EuMesh mesh(U, &rect);

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
    solver.prepare(mesh);
    double init_volume = volume(mesh);

    int n_step = 0;
    double end_time = 1.0;
    double curr_time = 0.0;
    double write_freq = end_time / 100;
    double write_next = 0.0;

    while (n_step < 10000) {
        if (curr_time >= write_next || curr_time >= end_time) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(8) << std::setprecision(3) << std::fixed
                      << curr_time << ";";

            double curr_volume = volume(mesh);
            std::cout << "\tLoss: " << std::setw(12) << std::setprecision(3) << std::scientific
                      << (curr_volume - init_volume) / (init_volume) << "\n";

            pvd.save(mesh.locals(), curr_time);

            AmrStorage crop = solver.body(mesh);
            pvd_body.save(crop, curr_time);

            AmrStorage exact = solver.exact(body, curr_time);
            pvd_exact.save(exact, curr_time);

            //AmrStorage scheme = solver.scheme(mesh);
            //pvd_scheme.save(scheme, curr_time);

            write_next += write_freq;

            if (curr_time >= end_time) {
                double vi = body.volume_inside(crop);
                std::cout << std::setprecision(3) << std::fixed;
                std::cout << "  Volume loss: " << 100 * (1.0 - vi / init_volume) << "%\n";
                pvd_body.save(crop, curr_time + 1.0e-13);
                break;
            }
        }

        // Определить шаг
        // double dt = solver.compute_dt(mesh); // для переменной скорости
        double dt = solver.get_dt();            // для постоянной скорости
        if (curr_time + dt > end_time) {
            dt = end_time - curr_time;
        }
        solver.set_dt(dt);

        // Шаг интегрирования
        if (splitting && solver.method() != Solver::Method::KABARE) {
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

void BodyLine::initial(Solver& solver, EuMesh& mesh, bool exact) const {
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

double BodyLine::volume_inside(AmrStorage &body) const {
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

void BodySquare::initial(Solver& solver, EuMesh& mesh, bool exact) const {
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

    if (solver.method() != Solver::Method::KABARE) {
        return;
    }

    for (int i = 0; i < mesh.nx(); ++i) {
        for (int j = 0; j < mesh.ny(); ++j) {
            auto cell = mesh(i, j);

            for (int k = 0; k < 4; k++) {
                cell(U).flow_u[k] = cell(U).u1;
                cell(U).flow_u_tmp[k] = cell(U).flow_u[k];
            }

            if (solver.velocity(cell.vs<+1, 0>()).x() > 0)  cell(U).flow_u[0] = cell(U).flow_u_tmp[0] =  cell(U).u1/solver.velocity(cell.vs<+1, 0>()).x();
            else if (solver.velocity(cell.vs<+1, 0>()).x() < 0)  cell(U).flow_u[0] = cell(U).flow_u_tmp[0] =  mesh(i+1,j).data(U).u1/solver.velocity(cell.vs<+1, 0>()).x();

            if (solver.velocity(cell.vs<-1, 0>()).x() > 0)  cell(U).flow_u[2] = cell(U).flow_u_tmp[2] =  mesh(i-1,j).data(U).u1/solver.velocity(cell.vs<-1, 0>()).x();
            else if (solver.velocity(cell.vs<+1, 0>()).x() < 0)  cell(U).flow_u[2] = cell(U).flow_u_tmp[2] =   cell(U).u1/solver.velocity(cell.vs<+1, 0>()).x();

            if (solver.velocity(cell.vs<0, 1>()).y() > 0)  cell(U).flow_u[1] = cell(U).flow_u_tmp[1] =  cell(U).u1/solver.velocity(cell.vs<0, 1>()).y();
            else if (solver.velocity(cell.vs<0, 1>()).y() < 0)  cell(U).flow_u[1] = cell(U).flow_u_tmp[1] =  mesh(i,j+1).data(U).u1/solver.velocity(cell.vs<0, 1>()).y();

            if (solver.velocity(cell.vs<0, -1>()).y() > 0)  cell(U).flow_u[3] = cell(U).flow_u_tmp[3] =  mesh(i,j-1).data(U).u1/solver.velocity(cell.vs<0, -1>()).y();
            else if (solver.velocity(cell.vs<0, 1>()).y() < 0)  cell(U).flow_u[3] = cell(U).flow_u_tmp[3] =  cell(U).u1/solver.velocity(cell.vs<0, 1>()).y();
        }
    }
}

double BodySquare::volume_inside(AmrStorage &body) const {
    double res = 0.0;
    for (auto& cell: body) {
        double vol = cell.approx_vol_fraction(inside);
        if (0.0 < vol && vol < 1.0) {
            vol = cell.volume_fraction(inside, 1000);
        }
        cell(U).u1 = vol;
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

void BodyDisk::initial(Solver& solver, EuMesh& mesh, bool exact) const {
    for (auto cell: mesh) {
        auto poly = cell.polygon();
        cell(U).u1 = poly.disk_clip_area(C, R) / cell.volume();
        cell(U).u2 = 0.0;
    }

    if (solver.method() != Solver::Method::KABARE) {
        return;
    }

    for (int i = 0; i < mesh.nx(); ++i) {
        for (int j = 0; j < mesh.ny(); ++j) {
            auto cell = mesh(i, j);

            for (int k = 0; k < 4; k++) {
                cell(U).flow_u[k] = cell(U).u1;
                cell(U).flow_u_tmp[k] = cell(U).flow_u[k];
            }

            if (solver.velocity(cell.vs<+1, 0>()).x() > 0)  cell(U).flow_u[0] = cell(U).flow_u_tmp[0] =  cell(U).u1;
            else if (solver.velocity(cell.vs<+1, 0>()).x() < 0)  cell(U).flow_u[0] = cell(U).flow_u_tmp[0] =  mesh(i+1,j).data(U).u1;

            if (solver.velocity(cell.vs<-1, 0>()).x() > 0)  cell(U).flow_u[2] = cell(U).flow_u_tmp[2] =  mesh(i-1,j).data(U).u1;
            else if (solver.velocity(cell.vs<+1, 0>()).x() < 0)  cell(U).flow_u[2] = cell(U).flow_u_tmp[2] =   cell(U).u1;

            if (solver.velocity(cell.vs<0, 1>()).y() > 0)  cell(U).flow_u[1] = cell(U).flow_u_tmp[1] =  cell(U).u1;
            else if (solver.velocity(cell.vs<0, 1>()).y() < 0)  cell(U).flow_u[1] = cell(U).flow_u_tmp[1] =  mesh(i,j+1).data(U).u1;

            if (solver.velocity(cell.vs<0, -1>()).y() > 0)  cell(U).flow_u[3] = cell(U).flow_u_tmp[3] =  mesh(i,j-1).data(U).u1;
            else if (solver.velocity(cell.vs<0, 1>()).y() < 0)  cell(U).flow_u[3] = cell(U).flow_u_tmp[3] =  cell(U).u1;
        }
    }
}

double BodyDisk::volume_inside(AmrStorage &body) const {
    double res = 0.0;
    for (auto& cell: body) {
        auto poly = cell.polygon();
        double vol = poly.disk_clip_area(C, R);
        cell(U).u1 = vol / cell.volume();
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

AmrStorage Solver::exact(Body& body, double curr_time) const {
    using zephyr::geom::Quad;
    using zephyr::geom::AmrCell;

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