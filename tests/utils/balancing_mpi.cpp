/// @brief Балансировка сетки

#include <iostream>
#include <zephyr/utils/mpi.h>

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

using zephyr::mesh::decomp::ORB;

struct _U_ {
    int rank;
    double load;
};

_U_ U;

double get_rank(const AmrStorage::Item &cell) { return cell(U).rank; }

double get_load(const AmrStorage::Item &cell) { return cell(U).load; }


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

std::vector<double> calc_loads(Mesh& mesh, int size) {
    std::vector<double> ws(size);
    for (auto cell: mesh) {
        ws[cell(U).rank] += cell(U).load;
    }
    return ws;
}

int main() {
    mpi::init();

    // Файл для записи
    PvdFile pvd("mesh", "output");

    pvd.variables += {"rank", get_rank};
    pvd.variables += {"load", get_load};

    // Сеточный генератор
    //Cuboid gen(0.0, 1.0, 0.0, 0.6, 0.0, 0.9);
    Rectangle gen(0.0, 1.0, 0.0, 0.6);
    gen.set_nx(200);

    // Bounding Box для сетки
    Box domain = gen.bbox();

    // Создаем сетку
    EuMesh mesh(U, gen);

    // Заполняем данные о нагрузке ячеек
    for (auto cell: mesh) {
        cell(U).load = foo(cell.center());
    }

    // Различные варианты инициализации ORB декомпозиции
    ORB orb(domain, "XY", mpi::size());
    mesh.set_decomposition(orb);
    //decomp::ORB orb(domain, "YX", 13);
    //decomp::ORB orb(domain, "YX", 13, 3);

    for (int step = 0; step < 1000; ++step) {
        // Вычисляем новый ранг ячеек
        for (auto &cell: mesh.locals()) {
            cell(U).rank = orb.rank(cell);
        }

        // Подсчитываем нагрузку каждого ранга
        auto ws = calc_loads(mesh, orb.size());

        if (step % 10 == 0) {
            std::cout << "Step:      " << step << "\n";
            std::cout << "Imbalance: " << ORB::imbalance(ws) << "\n";
            pvd.save(mesh, step);
        }

        // Балансировка декомпозиции
        orb.balancing(ws);
        mesh.redistribute();
    }

    mpi::finalize();
    return 0;
}