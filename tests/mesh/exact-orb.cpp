// Алгоритм точной ORB декомпозиции.

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

template <typename T>
double imbalance(const std::vector<T>& ws) {
    double max = *std::ranges::max_element(ws);
    double avg = std::accumulate(ws.begin(), ws.end(), 0.0) / static_cast<double>(ws.size());
    return  max / avg - 1.0;
}

int main() {
    mpi::handler init;

    // Сеточный генератор
    Cuboid gen(0.0, 1.0, 0.0, 0.6, 0.0, 0.5);
    gen.set_nx(70);
    //Rectangle gen(0.0, 1.5, 0.0, 1.0);
    //gen.set_nx(150);

    // Создаем сетку
    EuMesh mesh(gen);
    mesh.set_max_level(2);

    // Добавить переменные на сетку
    auto u = mesh.add<double>("u");

    // Bounding Box для сетки
    Box domain = gen.bbox();

    // Варианты инициализации ORB декомпозиции
    ORB orb(domain, "XYZ", mpi::size());
    //ORB orb(domain, "YX", 13);
    //ORB orb(domain, "YX", 13, 3);

    orb.use_exact(true);

    // Установить декомпозицию (+ делает redistribute)
    mesh.set_decomposition(orb);

    // Заполняем начальные данные и нагрузку
    auto init_cells = [u](EuCell& cell) {
        const Vector3d vc = {0.67, 0.43, 0.0};
        cell[u] = (cell.center() - vc).norm() < 0.23 ? 1.0 : 0.0;
    };

    auto set_flags = [u](EuCell& cell) {
        cell.set_flag(cell[u] > 0.5 ? 1 : -1);
    };

    for (int i = 0; i < mesh.max_level() + 2; ++i) {
        mesh.for_each(init_cells);
        mesh.for_each(set_flags);
        mesh.refine();
    }

    // Файл для записи
    PvdFile pvd("mesh", "output");
    pvd.variables = {"rank"};
    pvd.variables.append("u", u);

    pvd.save(mesh, 0.0);
    Stopwatch elapsed(true);
    mesh.balancing();
    elapsed.stop();
    mesh.redistribute(u);
    pvd.save(mesh, 1.0);

    auto ws = mpi::all_gather(static_cast<double>(mesh.n_cells()));
    mpi::cout << "Imbalance:  " << imbalance(ws) << "\n";

    mpi::cout << "\nBalancing time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    return 0;
}