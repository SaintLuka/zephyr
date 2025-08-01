/// @file advection.cpp
/// @brief Решение задачи переноса с неоднородной скоростью.
/// Используется простейшая схема первого порядка (upwind)

#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/pvd_file.h>

using zephyr::geom::Vector3d;
using zephyr::geom::Box;
using zephyr::geom::Boundary;
using zephyr::geom::generator::Rectangle;
using zephyr::mesh::Storable;
using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::io::PvdFile;

// Векторное поле скорости
Vector3d velocity(const Vector3d& c) {
    return { 1.0, 0.3 + 0.3*std::sin(4 * M_PI * c.x()), 0.0 };
}

int main() {
    // Геометрия области
    Rectangle rect(0.0, 1.0, 0.0, 0.6, true);
    rect.set_nx(200);
    rect.set_boundaries({
        .left   = Boundary::PERIODIC, .right = Boundary::PERIODIC,
        .bottom = Boundary::PERIODIC, .top   = Boundary::PERIODIC});

    // Создать сетку
    EuMesh mesh(rect);

    // Переменные для хранения на сетке
    auto u1 = mesh.add<double>("u1");
    auto u2 = mesh.add<double>("u2");

    // Файл для записи
    PvdFile pvd("mesh", "output");

    // Переменные для сохранения
    pvd.variables.append("u", u1);
    pvd.variables += {"vx", [](EuCell cell) -> double { return velocity(cell.center()).x(); } };
    pvd.variables += {"vy", [](EuCell cell) -> double { return velocity(cell.center()).y(); } };

    // Заполняем начальные данные
    Box box = mesh.bbox();
    Vector3d vc = box.center();
    double D = 0.1 * box.diameter();
    for (auto cell: mesh) {
        cell(u1) = (cell.center() - vc).norm() < D ? 1.0 : 0.0;
        cell(u2) = 0.0;
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
            pvd.save(mesh, curr_time);
            next_write += 0.02;
        }

        // Определяем dt
        double dt = std::numeric_limits<double>::max();
        for (auto cell: mesh) {
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
            double zc = cell(u1);

            double fluxes = 0.0;
            for (auto& face: cell.faces()) {
                double zn = face.neib(u1);

                double af = velocity(face.center()).dot(face.normal());
                double a_p = std::max(af, 0.0);
                double a_m = std::min(af, 0.0);

                fluxes += (a_p * zc + a_m * zn) * face.area();
            }

            cell(u2) = zc - dt * fluxes / cell.volume();
        }

        // Обновляем слои
        for (auto cell: mesh) {
            cell(u1) = cell(u2);
            cell(u2) = 0.0;
        }

        n_step += 1;
        curr_time += dt;
    }

    return 0;
}