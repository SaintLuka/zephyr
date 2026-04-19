// Алгоритм точной ORB декомпозиции.

#include <iostream>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/mesh/decomp/ORB.h>

#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/rectangle.h>

#include "zephyr/geom/generator/array2d.h"
#include "zephyr/utils/pyplot.h"

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

void test_points_single() {
    int size = 9; // число блоков декомпозиции

    std::vector<Vector3d> points(53 * size);

    for (auto& p : points) {
        p.x() = rand() / double(RAND_MAX);
        p.y() = rand() / double(RAND_MAX);
    }

    Box domain = Box::Empty(2);
    for (auto& p: points) {
        domain.capture(p);
    }


    ORB orb(domain, "X", size);

    std::vector<int> ranks(points.size());
    std::vector<double> f_ranks(points.size());
    for (int i = 0; i < points.size(); ++i) {
        ranks[i] = orb.rank(points[i]);
        f_ranks[i] = ranks[i];
    }

    std::vector<int> block_count(size, 0);
    for (auto r: ranks) {
        ++block_count[r];
    }

    std::cout << "Counts: ";
    for (auto r: block_count) {
        std::cout << r << ", ";
    }
    std::cout << "\nImbalance: " << imbalance(block_count) << "\n";


    orb.exact_balancing(points);


    for (int i = 0; i < points.size(); ++i) {
        ranks[i] = orb.rank(points[i]);
        f_ranks[i] = ranks[i];
    }

    std::ranges::fill(block_count, 0);
    for (auto r: ranks) {
        ++block_count[r];
    }

    std::cout << "Counts: ";
    for (auto r: block_count) {
        std::cout << r << ", ";
    }
    std::cout << "\nImbalance: " << imbalance(block_count) << "\n";



    zephyr::utils::pyplot plt;
    plt.set_aspect_equal();
    std::vector<double> xs(points.size());
    std::vector<double> ys(points.size());
    for (int i = 0; i < points.size(); ++i) {
        xs[i] = points[i].x();
        ys[i] = points[i].y();
    }
    plt.scatter(xs, ys, {.c=f_ranks, .cmap="Paired"});

    auto lines = orb.blocks().lines();
    for (const auto& line: lines) {
        xs.resize(line.size());
        ys.resize(line.size());
        for (int i = 0; i < line.size(); ++i) {
            xs[i] = line[i].x();
            ys[i] = line[i].y();
        }
        plt.plot(xs, ys, {.color="black"});
    }

    plt.tight_layout();
    plt.show();



}

void test_points_mpi() {
    throw std::runtime_error("Not implemented");
}

void test_points() {
    if (mpi::single()) {
        test_points_single();
    }
    else {
        test_points_mpi();
    }
}

void test_grid() {
    // Сеточный генератор
    //Cuboid gen(0.0, 1.0, 0.0, 0.6, 0.0, 0.9);
    //gen.set_nx(50);
    Rectangle gen(0.0, 1.5, 0.0, 1.0);
    gen.set_nx(150);

    // Создаем сетку
    EuMesh mesh(gen);
    mesh.set_max_level(3);

    // Добавить переменные на сетку
    auto u = mesh.add<double>("u");

    // Bounding Box для сетки
    Box domain = gen.bbox();

    // Варианты инициализации ORB декомпозиции
    ORB orb(domain, "X", mpi::size());
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
    mesh.redistribute(u);
    elapsed.stop();
    pvd.save(mesh, 1.0);

    auto ws = mpi::all_gather(static_cast<double>(mesh.n_cells()));
    mpi::cout << "Imbalance:  " << mesh.decomp().imbalance(ws) << "\n";

    mpi::cout << "\nElapsed time:   " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";
}

int main(int argc, char** argv) {
    mpi::handler init(argc, argv);

    test_points();

    return 0;
}