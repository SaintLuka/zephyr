/// @file transfer_3D.cpp
/// @brief Решение задачи переноса со сложными решателями.
#include <iostream>
#include <iomanip>

#include <zephyr/io/pvd_file.h>
#include <zephyr/geom/geom.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/surface/solid_body_3d.h>
#include <zephyr/math/solver/transfer.h>

using namespace zephyr::geom;
using namespace zephyr::mesh;

using zephyr::io::PvdFile;
using generator::Cuboid;

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
            Vector3d V0 = {0.7, -0.35, 0.3};
            return V0;
        }
        if (test == Test::Rotation) {
            Vector3d center = {0.5, 0.5, 0.3}; // Центр вращения
            Vector3d axis = {-0.3, 0.1, 0.7};  // Ось вращения
            axis = M_PI * axis.normalized();
            return axis.cross(p - center);
        }
        return Vector3d::Zero();
    }

    // Сетка с точным решением от времени
    Mesh exact(SolidBody3D& body, double curr_time) const {
        // Точное решение
        if (test == Test::Translation) {
            Vector3d V0 = {0.7, -0.35, 0.3};
            body.move(curr_time * V0);
        }
        else if (test == Test::Rotation) {
            Vector3d center = {0.5, 0.5, 0.3}; // Центр вращения
            Vector3d axis = {-0.3, 0.1, 0.7};  // Ось вращения
            body.rotation_relative(axis, M_PI * curr_time, center);
        }

        auto triangles = body.triangulation(400);
        int n_triangles = triangles.size();

        Mesh cells = Mesh::PolySet(2);
        cells.locals().reserve(n_triangles, 3 * n_triangles, 3 * n_triangles);
        for (const auto& tri: triangles) {
            Polygon poly(tri);
            cells.push_back(poly);
        }
        return cells;
    }
};

static Solver::State data;

// Объем тела
double volume(Mesh& cells, Storable<double> u1) {
    double sum = 0.0;
    for (auto cell: cells) {
        sum += cell.volume() * cell[u1];
    }
    return sum;
}

// Какой объем сетки отсекается телом
double volume_inside(const SolidBody3D& body, Mesh &cells) {
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
    pvd_body.options.polyhedral = true;

    // Геометрия области
    Cuboid gen(0.0, 1.0, 0.0, 0.7, 0.0, 0.6);
    gen.set_nx(50);
    gen.set_boundaries({
        .left   = Boundary::ZOE, .right = Boundary::ZOE,
        .bottom = Boundary::ZOE, .top   = Boundary::ZOE,
        .back   = Boundary::ZOE, .front = Boundary::ZOE});

    // Создать сетку
    Mesh mesh(gen);
    mesh.set_max_level(2);

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
    BodyBall body(0.1, {0.15, 0.5, 0.15});
    //BodyCube body(0.15, {0.15, 0.5, 0.15});

    // Выбора теста
    solver.test = Solver::Test::Rotation;

    // Переменные для сохранения
    pvd.variables = {"level"};
    pvd.variables.add_cell_data("u", data.u1);
    pvd.variables.add_cell_data("u2", data.u2);
    pvd.variables.add_cell_data("n", data.n);
    pvd.variables.add_cell_data("p", data.p);
    pvd.variables.add_cell_data("grad", data.grad);
    pvd.variables += {"over", [](Cell& cell) -> double {
        double u = cell[data.u1];
        return u < 0.0 ? u : (u <= 1.0 ? NAN : u - 1.0);
    }};
    pvd.variables += {"close", [](Cell& cell) -> double {
        double u = cell[data.u1];
        return std::abs(u < 0.5 ? u : 1.0 - u);
    }};
    pvd.variables += {"min_val", [](Cell& cell) -> double {
        double min_val = cell[data.u1];
        for (auto face: cell.faces()) {
            if (face.is_boundary()) continue;
            min_val = std::min(min_val, face.neib(data.u1));
        }
        return min_val;
    }};
    pvd.variables += {"max_val", [](Cell& cell) -> double {
        double max_val = cell[data.u1];
        for (auto face: cell.faces()) {
            if (face.is_boundary()) continue;
            max_val = std::max(max_val, face.neib(data.u1));
        }
        return max_val;
    }};

    // Начальные условия
    auto initialize = [&body](Cell& cell) {
        cell[data.u1] = body.volume_fraction(cell, 1.0e-4);
        cell[data.u2] = 0.0;
    };

    // Выполняет инициализацию при адаптации
    Distributor init_distr = Distributor::initializer(initialize);

    std::cout << "Initialization\n";
    mesh.for_each(initialize);
    mesh.set_distributor(init_distr);
    for (int lvl = 0; lvl <= mesh.max_level(); ++lvl) {
        std::cout << "    Level " << lvl << " / " << mesh.max_level() << "\n";
        solver.set_flags(mesh);
        mesh.refine();
    }
    mesh.set_distributor(solver.distributor());

    solver.update_interface(mesh);
    double init_volume = volume(mesh, data.u1);

    int n_step = 0;
    double end_time = 1.0;
    double curr_time = 0.0;
    double write_freq = end_time / 50;
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

            solver.update_interface(mesh);
            Mesh crop = solver.body(mesh);
            pvd_body.save(crop, curr_time);

            auto curr_body = body;
            Mesh exact = solver.exact(curr_body, curr_time);
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
            solver.update(mesh, Direction::Z);
        }
        else {
            solver.update(mesh, Direction::ANY);
        }

        // Адаптация
        solver.set_flags(mesh);
        mesh.refine();

        n_step += 1;
        curr_time += solver.get_dt();
    }

    // Финальные записи
    pvd.save(mesh, curr_time);

    solver.update_interface(mesh);
    Mesh crop = solver.body(mesh);
    pvd_body.save(crop, curr_time);

    Mesh exact = solver.exact(body, curr_time);
    pvd_exact.save(exact, curr_time);

    double vi = volume_inside(body, crop);
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\n  Final mismatch: " << 100 * (1.0 - vi / init_volume) << "%\n";
    pvd_body.save(crop, curr_time + 1.0e-13);

    return 0;
}
