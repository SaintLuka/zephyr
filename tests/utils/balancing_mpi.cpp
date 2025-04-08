/// @brief Балансировка сетки на mpi
//
// TODO: удалить, плохая копия distributed-conv

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
    int rank;
    double load;
    double u1, u2;
};

_U_ U;

// Векторное поле скорости
Vector3d velocity(const Vector3d& c) {
    return { 1.0, 0.3 + 0.3*std::sin(4 * M_PI * c.x()), 0.0 };
}

double get_rank(const AmrStorage::Item &cell) { return cell(U).rank; }

double get_load(const AmrStorage::Item &cell) { return cell(U).load; }

// Переменные для сохранения
double get_u(AmrStorage::Item& cell)  {
    return cell(U).u1;
}

double get_vx(AmrStorage::Item& cell) {
    return velocity(cell.center).x();
}

double get_vy(AmrStorage::Item& cell) {
    return velocity(cell.center).y();
}

double foo(const Vector3d& v) {
    Vector3d c1 = {0.2, 0.4, 0.0};
    Vector3d c2 = {0.6, 0.1, 0.0};
    Vector3d c3 = {0.8, 0.4, 0.0};
    double A1 = 200.0;
    double A2 = 100.0;
    double A3 = 50.0;
    double s1 = 0.15;
    double s2 = 0.2;
    double s3 = 0.1;
    return A1/s1 * std::exp(-(v - c1).squaredNorm() / (s1 * s1)) +
           A2/s2 * std::exp(-(v - c2).squaredNorm() / (s2 * s2)) +
           A3/s3 * std::exp(-(v - c3).squaredNorm() / (s3 * s3));
}

double calc_loads(Mesh& mesh, int size) {
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

    pvd.variables += {"rank", get_rank};
    pvd.variables += {"load", get_load};
    pvd.variables += {"u",  get_u};
    pvd.variables += {"vx", get_vx};
    pvd.variables += {"vy", get_vy};

    // Сеточный генератор
    //Cuboid gen(0.0, 1.0, 0.0, 0.6, 0.0, 0.9);
    Rectangle gen(0.0, 1.0, 0.0, 0.6);
    gen.set_nx(200);

    // Bounding Box для сетки
    Box domain = gen.bbox();

    // Создаем сетку
    EuMesh mesh(gen, U);

    // Заполняем данные о нагрузке ячеек
    Vector3d vc = domain.center();
    double D = 0.1 * domain.diameter();
    for (auto cell: mesh) {
        cell(U).u1 = (cell.center() - vc).norm() < D ? 1.0 : 0.0;
        cell(U).u2 = 0.0;
        cell(U).load = std::sqrt(cell(U).u1 * cell(U).u1 + cell(U).u2 * cell(U).u2);
    }

    // Различные варианты инициализации ORB декомпозиции
    ORB orb(domain, "XY", mpi::size());
    mesh.set_decomposition(orb);
    //decomp::ORB orb(domain, "YX", 13);
    //decomp::ORB orb(domain, "YX", 13, 3);

    // Число Куранта
    double CFL = 0.5;

    int n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    Stopwatch elapsed(true);
    for (int step = 0; step < 1000; ++step) {
        // Балансировка декомпозиции
        {
            double load = calc_loads(mesh, orb.size());
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
            cell(U).load = std::sqrt(cell(U).u1 * cell(U).u1 + cell(U).u2 * cell(U).u2);
        }

        n_step += 1;
        curr_time += dt;
    }
    elapsed.stop();

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    return 0;
}