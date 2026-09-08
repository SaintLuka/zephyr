/// @file transfer_2D.cpp
/// @brief Решение задачи переноса со сложными решателями.
#include <iostream>
#include <iomanip>

#include <zephyr/io/pvd_file.h>
#include <zephyr/geom/geom.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/surface/solid_body_2d.h>
#include <zephyr/math/solver/transfer.h>

using namespace zephyr::geom;
using namespace zephyr::mesh;

using zephyr::io::PvdFile;
using generator::Rectangle;

// Наследуем собственный решатель от Transfer, теперь переопределив
// поле скорости можно решать произвольные задачи на перенос.
class Solver : public zephyr::math::Transfer {
public:
    enum class Test {
        Translation,  // Поступательное движение
        Rotation,     // Вращательное движение
    };

    Test test = Test::Translation;

    // Скорость переноса в соответствии с Solver::test
    Vector3d velocity(const Vector3d& p) const override {
        if (test == Test::Translation) {
            Vector3d V0 = {0.7, -0.35, 0.0};
            return V0;
        }
        if (test == Test::Rotation) {
            Vector3d center = {0.5, 0.5, 0.0};  // Центр вращения
            Vector3d omega = {0.0, 0.0, M_PI};  // Угловая частота
            return omega.cross(p - center);
        }
        return Vector3d::Zero();
    }

    // Сетка с точным решением от времени
    EuMesh exact(SolidBody2D& body, double curr_time) const {
        // Точное решение
        if (test == Test::Translation) {
            Vector3d V0 = {0.7, -0.35, 0.0};
            body.move(curr_time * V0);
        }
        else if (test == Test::Rotation) {
            double omega = M_PI;  // Угловая частота
            Vector3d center = {0.5, 0.5, 0.0};  // Центр вращения
            body.rotation_relative(omega * curr_time, center);
        }

        auto vs = body.outline(100);
        EuMesh cells = EuMesh::PolySet(2);
        for (size_t i = 0; i < vs.size(); ++i) {
            size_t j = (i + 1) % vs.size();
            Line line = {vs[i], vs[j]};
            cells.push_back(line);
        }
        return cells;
    }
};

static Solver::State data;

// Объем тела
double volume(EuMesh& cells, Storable<double> u1) {
    double sum = 0.0;
    for (auto cell: cells) {
        sum += cell.volume() * cell[u1];
    }
    return sum;
}

// Какой объем сетки отсекается телом
double volume_inside(const SolidBody2D& body, EuMesh &cells) {
    double sum = 0.0;
    for (auto cell: cells) {
        sum += body.volume_inside(cell, 1.0e-3);
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

    // Создать сетку
    EuMesh mesh(rect);

    // Создать решатель
    Solver solver;

    // Добавить типы
    data = solver.add_types(mesh);

    // Настройки решателя
    solver.set_CFL(0.5);
    solver.set_dim(mesh.dim());
    solver.set_plic_type(Plic::CSIR);

    // CRP_V3 сейчас не как на картинке, остальное всё повторяется
    solver.set_method(Solver::Method::CRP_N2);

    // Выбор фигуры
    //BodyStrip body(0.1, {0.15, 0.5, 0.0});
    //BodyDisk body(0.1, {0.15, 0.5, 0.0});
    BodySquare body(0.2, {0.15, 0.5, 0.0});

    // Выбор теста
    solver.test = Solver::Test::Rotation;

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables.append("u", data.u1);
    pvd.variables.append("u2", data.u2);
    pvd.variables.append("n", data.n);
    pvd.variables.append("p", data.p);
    pvd.variables.append("grad", data.grad);
    pvd.variables += {"over", [](EuCell& cell) -> double {
        double u = cell[data.u1];
        return u < 0.0 ? u : (u <= 1.0 ? NAN : u - 1.0);
    }};
    pvd.variables += {"close", [](EuCell& cell) -> double {
        double u = cell[data.u1];
        return std::abs(u < 0.5 ? u : 1.0 - u);
    }};

    // Начальные условия
    for (auto cell: mesh) {
        cell[data.u1] = body.volume_fraction(cell, 1.0e-4);
        cell[data.u2] = 0.0;
    }

    solver.update_interface(mesh);
    double init_volume = volume(mesh, data.u1);

    int n_step = 0;
    double end_time = 1.0;
    double curr_time = 0.0;
    double write_freq = end_time / 100;
    double write_next = 0.0;

    // Расщепление по направлениям
    bool splitting = true;

    while (curr_time < end_time) {
        if (curr_time >= write_next || curr_time >= end_time) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(8) << std::setprecision(3) << std::fixed
                      << curr_time << ";";

            double curr_volume = volume(mesh, data.u1);
            std::cout << "\tLoss: " << std::setw(12) << std::setprecision(3) << std::scientific
                      << (curr_volume - init_volume) / (init_volume) << "\n";

            pvd.save(mesh, curr_time);

            EuMesh crop = solver.body(mesh);
            pvd_body.save(crop, curr_time);

            auto curr_body = body;
            EuMesh exact = solver.exact(curr_body, curr_time);
            pvd_exact.save(exact, curr_time);

            write_next += write_freq;
        }

        // Точное завершение в end_time
        solver.set_max_dt(end_time - curr_time);

        // Определить шаг
        double dt = solver.compute_dt(mesh);
        solver.set_dt(dt);

        // Шаг интегрирования
        if (splitting) {
            solver.update(mesh, Direction::X);
            solver.update(mesh, Direction::Y);
        }
        else {
            solver.update(mesh, Direction::ANY);
        }

        n_step += 1;
        curr_time += solver.get_dt();
    }

    // Финальные записи
    pvd.save(mesh, curr_time);

    EuMesh crop = solver.body(mesh);
    pvd_body.save(crop, curr_time);

    EuMesh exact = solver.exact(body, curr_time);
    pvd_exact.save(exact, curr_time);

    double vi = volume_inside(body, crop);
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\n  Final mismatch: " << 100 * (1.0 - vi / init_volume) << "%\n";
    pvd_body.save(crop, curr_time + 1.0e-13);

    return 0;
}
