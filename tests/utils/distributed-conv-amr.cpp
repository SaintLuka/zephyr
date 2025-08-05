// @brief Динамическая ORB декомпозиция сетки совместно с адаптацией.
// На сетке решается уравнение переноса для проверки связности сетки.

#include <iostream>
#include <iomanip>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/mesh/decomp/ORB.h>

#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/rectangle.h>

using namespace zephyr::mesh;
using namespace zephyr::geom;

using zephyr::utils::mpi;
using zephyr::utils::Stopwatch;
using zephyr::io::PvdFile;
using zephyr::geom::generator::Cuboid;
using zephyr::geom::generator::Rectangle;

using zephyr::mesh::decomp::ORB;


// Векторное поле скорости
Vector3d velocity(const Vector3d& c) {
    return { 1.0, 0.3 + 0.3*std::sin(4 * M_PI * c.x()), 0.3 };
}

// Шар в левом нижнем углу области
void set_initials(EuMesh& mesh, const Box& domain,
        Storable<double> u1, Storable<double> u2) {
    Vector3d vc = domain.vmin + 0.2 * domain.size();
    double D = 0.1 * domain.diameter();
    for (auto cell: mesh) {
        cell(u1) = (cell.center() - vc).norm() < D ? 1.0 : 0.0;
        cell(u2) = 0.0;
    }
}

// Адаптация больших перепадов концентрации
void set_flags(EuMesh& mesh, Storable<double> var) {
    for (auto cell: mesh) {
        cell.set_flag(-1);

        double uc = cell(var);
        for (auto face: cell.faces()) {
            double un = face.neib(var);

            if (std::abs(uc - un) > 0.01) {
                cell.set_flag(1);
                break;
            }
        }
    }
}

// Распределитель данных при адаптации
Distributor get_distributor(Storable<double> var) {
    Distributor distr = Distributor::simple();
    distr.merge = [var](Children &children, EuCell &parent) {
        double sum = 0.0;
        for (auto child: children) {
            sum += child(var) * child.volume();
        }
        parent(var) = sum / parent.volume();
    };
    return distr;
}

int main() {
    mpi::handler init;

    // Сеточный генератор
    //Cuboid gen(0.0, 1.0, 0.0, 0.6, 0.0, 0.9);
    Rectangle gen(0.0, 1.0, 0.0, 1.0);
    gen.set_nx(123);

    // Создаем сетку
    EuMesh mesh(gen);

    // Добавить данные на сетку
    auto [u1, u2] = mesh.add<double>("u1", "u2");

    mesh.set_max_level(3);
    mesh.set_distributor(get_distributor(u1));

    // Добавляем декомпозицию
    mesh.set_decomposition("XY");

    // Bounding Box для сетки
    Box domain = gen.bbox();

    // Начальные данные
    for (int i = 0; i < mesh.max_level(); ++i) {
        set_initials(mesh, domain, u1, u2);
        mesh.sync(u1);
        set_flags(mesh, u1);
        mesh.refine();
    }
    set_initials(mesh, domain, u1, u2);

    // Файл для записи
    PvdFile pvd("mesh", "output");

    pvd.variables = {"rank", "level"};
    pvd.variables.append("u", u1);

    // Число Куранта
    double CFL = 0.5;

    int n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    Stopwatch elapsed(true);
    while (curr_time <= 0.8) {
        // Балансировка декомпозиции
        if (n_step % 10 == 0) {
            mesh.balancing(mesh.n_cells());
            mesh.redistribute(u1);
        }

        if (curr_time >= next_write) {
            mpi::cout << "\tШаг: " << std::setw(6) << n_step << ";"
                      << "\tВремя: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);
            next_write += 0.02;
        }

        // Определяем dt
        double dt = std::numeric_limits<double>::max();
        for (auto& cell: mesh) {
            double dx = cell.incircle_diameter();
            dt = std::min(dt, dx / velocity(cell.center()).norm());
        }
        dt = mpi::min(CFL * dt);

        // Отправить/получить основной слой
        mesh.sync(u1);

        // Расчет потоков по схеме upwind
        for (auto cell: mesh) {
            double fluxes = 0.0;
            for (auto& face: cell.faces()) {
                auto neib = face.neib();

                double af = velocity(face.center()).dot(face.normal());
                double a_p = std::max(af, 0.0);
                double a_m = std::min(af, 0.0);

                fluxes += (a_p * cell(u1) + a_m * neib(u1)) * face.area();
            }

            cell(u2) = cell(u1) - dt * fluxes / cell.volume();
        }

        // Обновляем слои
        for (auto& cell: mesh) {
            cell(u1) = cell(u2);
            cell(u2) = 0.0;
        }

        // Отправить/получить основной слой
        mesh.sync(u1);
        set_flags(mesh, u1);
        mesh.refine();

        n_step += 1;
        curr_time += dt;
    }
    elapsed.stop();

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    return 0;
}