// @brief Динамическая ORB декомпозиция без адаптации.
// Можно проверить динамическую балансировку на различных сетках,
// включая неструктурированные сетки из многоугольников.

#include <iostream>

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
using zephyr::io::PvdFile;
using zephyr::geom::generator::Cuboid;
using zephyr::geom::generator::Rectangle;
using zephyr::utils::Stopwatch;

using zephyr::mesh::decomp::ORB;

// Векторное поле скорости
Vector3d velocity(const Vector3d& c) {
    return { 1.0, 0.5 + 0.3*std::sin(4 * M_PI * c.x()), 0.0 };
}

// Просуммировать фиктивную нагрузку
double calc_loads(EuMesh& mesh, Storable<double> load) {
    double full = 0;
    for (auto cell: mesh) {
        full += cell(load);
    }
    return full;
}

int main() {
    mpi::handler init;

    // Сеточный генератор
    //Cuboid gen(0.0, 1.0, 0.0, 0.6, 0.0, 0.9);
    //gen.set_nx(50);
    Rectangle gen(0.0, 1.0, 0.0, 1.0, true);
    gen.set_nx(400);

    // Создаем сетку
    EuMesh mesh(gen);

    // Добавить переменные на сетку
    auto [u1, u2, load] = mesh.add<double>("u1", "u2", "load");

    // Bounding Box для сетки
    Box domain = gen.bbox();

    // Варианты инициализации ORB декомпозиции
    ORB orb(domain, "XY", mpi::size());
    //ORB orb(domain, "YX", 13);
    //ORB orb(domain, "YX", 13, 3);

    // Установить декомпозицию (+ делает redistribute)
    mesh.set_decomposition(orb);

    // Заполняем начальные данные и нагрузку
    Vector3d vc = domain.vmin + 0.2 * domain.size();
    double D = 0.1 * domain.diameter();
    for (auto cell: mesh) {
        cell(u1) = (cell.center() - vc).norm() < D ? 1.0 : 0.0;
        cell(u2) = 0.0;
        cell(load) = cell(u1);
    }

    // Файл для записи
    PvdFile pvd("mesh", "output");

    pvd.variables = {"rank", "index"};
    pvd.variables.append("u", u1);
    pvd.variables.append("load", load);

    // Число Куранта
    double CFL = 0.5;

    Stopwatch elapsed(true);
    for (int step = 0; step < 200; ++step) {
        // Балансировка декомпозиции
        {
            double full = calc_loads(mesh, load);
            mesh.balancing(full);
            mesh.redistribute(u1);
        }

        // Подсчитываем нагрузку каждого ранга
        if (step % 10 == 0) {
            mpi::cout << "Step:      " << step << "\n";
            pvd.save(mesh, step);
        }

        double dt = std::numeric_limits<double>::max();
        for (auto& cell: mesh) {
            double dx = cell.incircle_diameter();
            dt = std::min(dt, dx / velocity(cell.center()).norm());
        }
        dt = mpi::min(CFL * dt);

        // Отправить/получить основной слой
        mesh.sync(u1);

        // Расчет по схеме upwind
        for (auto& cell: mesh) {
            double zc = cell(u1);

            double fluxes = 0.0;
            for (auto& face: cell.faces()) {
                double zn = face.neib(u1);

                double af = velocity(face.center()).dot(face.normal());
                double a_p = std::max(af, 0.0);
                double a_m = std::min(af, 0.0);

                fluxes += (a_p * zc + a_m * zn) * face.area();
            }

            cell(u2) = cell(u1) - dt * fluxes / cell.volume();
        }

        // Обновляем слои
        for (auto& cell: mesh) {
            cell(u1) = cell(u2);
            cell(u2) = 0.0;
            cell(load) = cell(u1);
        }
    }
    elapsed.stop();

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    return 0;
}