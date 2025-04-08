/// @file Динамическая ORB декомпозиция без адаптации.
/// Можно проверить динамическую балансировку на различных сетках,
/// включая неструктурированные сетки из многоугольников.

#include <iostream>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/mesh/decomp/ORB.h>

#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/rectangle.h>

using namespace zephyr::mesh;
using zephyr::utils::mpi;
using zephyr::io::PvdFile;
using zephyr::geom::generator::Cuboid;
using zephyr::geom::generator::Rectangle;
using zephyr::utils::Stopwatch;

using zephyr::mesh::decomp::ORB;

struct _U_ {
    double load;
    double u1, u2;
};

_U_ U;

// Векторное поле скорости
Vector3d velocity(const Vector3d& c) {
    return { 1.0, 0.3 + 0.3*std::sin(4 * M_PI * c.x()), 0.0 };
}

// Переменные для сохранения
double get_u(AmrStorage::Item& cell)  { return cell(U).u1; }

double get_load(const AmrStorage::Item &cell) { return cell(U).load; }

// Просуммировать фиктивную нагрузку
double calc_loads(Mesh& mesh) {
    double load = 0;
    for (auto cell: mesh) {
        load += cell(U).load;
    }
    return load;
}

int main() {
    mpi::handler init;

    // Файл для записи
    PvdFile pvd("mesh", "output");

    pvd.variables = {"rank"};
    pvd.variables += {"u",  get_u};
    pvd.variables += {"load", get_load};

    // Сеточный генератор
    //Cuboid gen(0.0, 1.0, 0.0, 0.6, 0.0, 0.9);
    //gen.set_nx(50);
    Rectangle gen(0.0, 1.0, 0.0, 1.0, true);
    gen.set_nx(250);

    // Создаем сетку
    EuMesh mesh(gen, U);

    // Bounding Box для сетки
    Box domain = gen.bbox();

    // Заполняем данные о нагрузке ячеек
    Vector3d vc = domain.vmin + 0.2 * domain.size();
    double D = 0.1 * domain.diameter();
    for (auto cell: mesh) {
        cell(U).u1 = (cell.center() - vc).norm() < D ? 1.0 : 0.0;
        cell(U).u2 = 0.0;
        cell(U).load = cell(U).u1;
    }

    // Различные варианты инициализации ORB декомпозиции
    ORB orb(domain, "XY", mpi::size());
    //ORB orb(domain, "YX", 13);
    //ORB orb(domain, "YX", 13, 3);

    mesh.set_decomposition(orb);

    // Число Куранта
    double CFL = 0.5;

    int n_step = 0;

    Stopwatch elapsed(true);
    for (int step = 0; step < 500; ++step) {
        // Балансировка декомпозиции
        {
            double load = calc_loads(mesh);
            mesh.balancing(load);
            mesh.redistribute();
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

        // Расчет по схеме upwind
        for (auto& cell: mesh) {
            auto& zc = cell(U);

            double fluxes = 0.0;
            for (auto& face: cell.faces()) {
                const auto& zn = face.neib(U);

                double af = velocity(face.center()).dot(face.normal());
                double a_p = std::max(af, 0.0);
                double a_m = std::min(af, 0.0);

                fluxes += (a_p * zc.u1 + a_m * zn.u1) * face.area();
            }

            zc.u2 = zc.u1 - dt * fluxes / cell.volume();
        }

        // Обновляем слои
        for (auto& cell: mesh) {
            cell(U).u1 = cell(U).u2;
            cell(U).u2 = 0.0;
            cell(U).load = cell(U).u1;
        }

        n_step += 1;
    }
    elapsed.stop();

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    return 0;
}