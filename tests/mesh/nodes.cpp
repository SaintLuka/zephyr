#include <iostream>
#include <iomanip>

#include <zephyr/io/pvd_file.h>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/threads.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/geom/boundary.h>
#include <zephyr/geom/grid.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/collection/plane_with_hole.h>

#include <zephyr/mesh/euler/eu_prim.h>
#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/mesh/decomp/ORB.h>
#include <zephyr/mesh/decomp/rwalk.h>

using namespace zephyr::geom;
using namespace zephyr::mesh;
using namespace zephyr::io;
using namespace zephyr::geom::generator;

using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::mesh::decomp::ORB;
using zephyr::mesh::decomp::RWalk;
using zephyr::utils::mpi;
using zephyr::utils::threads;

auto WALL = Boundary::WALL;
auto ZOE = Boundary::ZOE;


// Простой квадрат с декартовой сеткой
EuMesh test1() {
    Rectangle gen(-1.0, 1.0, -1.0, 1.0);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_nx(20);
    gen.set_adaptive(false);
    return EuMesh(gen);
}

// Простой квадрат с адаптивной декартовой сеткой
EuMesh test2() {
    Rectangle gen(-1.0, 1.0, -1.0, 1.0);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_nx(20);
    return EuMesh(gen);
}

// Ячейки Вороного в прямоугольнике
EuMesh test3() {
    Rectangle gen(-1.0, 1.0, -1.0, 1.0, true);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_nx(20);
    return EuMesh(gen);
}

// Декартова сетка в кубе
EuMesh test4() {
    Cuboid gen(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL, .back=WALL, .front=WALL});
    gen.set_nx(20);
    gen.set_adaptive(false);
    return EuMesh(gen);
}

// Декартова адаптивная сетка в кубе
EuMesh test5() {
    Cuboid gen(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL, .back=WALL, .front=WALL});
    gen.set_nx(20);
    return EuMesh(gen);
}

// Extrude сетки из многоугольников
EuMesh test6() {
    Rectangle gen(-1.0, 1.0, -1.0, 1.0, true);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_nx(20);
    Grid grid = gen.make();
    grid.extrude(Vector3d{0.0, 0.0, 1.0}, 12, WALL, WALL);
    return EuMesh(std::move(grid));
}

// Простая BlockStructured сетка
EuMesh test7() {
    collection::PlaneWithHole gen(0.0, 2.0, 0.0, 1.0, 0.6, 0.2, 0.1);
    gen.set_boundaries({.left=ZOE, .right=ZOE, .bottom=WALL, .top=WALL});
    gen.set_ny(40);
    return EuMesh(gen);
}

// BlockStructured с extrude
EuMesh test8() {
    collection::PlaneWithHole gen = collection::PlaneWithHole(0.0, 2.0, 0.0, 2.0, 1.0, 1.0, 0.3);
    gen.set_boundaries({.left=ZOE, .right=ZOE, .bottom=WALL, .top=WALL});
    gen.set_ny(40);
    Grid grid = gen.make();
    grid.extrude(Vector3d::UnitZ(), 20, ZOE, ZOE);
    return EuMesh(std::move(grid));
}

// BlockStructured, затем make_amr
EuMesh test9() {
    collection::PlaneWithHole gen = collection::PlaneWithHole(0.0, 2.0, 0.0, 2.0, 1.0, 1.0, 0.3);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_ny(40);
    Grid grid = gen.make();
    grid.make_amr();
    return EuMesh(std::move(grid));
}

// BlockStructured, затем extrude и make_amr
EuMesh test10() {
    collection::PlaneWithHole gen = collection::PlaneWithHole(0.0, 2.0, 0.0, 2.0, 1.0, 1.0, 0.3);
    gen.set_boundaries({.left=WALL, .right=WALL, .bottom=WALL, .top=WALL});
    gen.set_ny(40);
    Grid grid = gen.make();
    grid.extrude(Vector3d::UnitZ()/5, 4, ZOE, ZOE);
    grid.make_amr();
    return EuMesh(std::move(grid));
}

// План.
// make_unique_nodes()
// check_base() for all
// redistribute
// check_base() for all

void check_mesh(const EuMesh& mesh) {
    mpi::for_each([&]() {
        int res = mesh.check_base();
        if (res < 0) {
            std::cout << "  rank " << mpi::rank() << ": check mesh failed!\n";
        }
        else {
            std::cout << "  rank " << mpi::rank() << ": mesh is fine!\n";
        }
    });
}

void save_markers(const AmrNodes& nodes, std::string filename) {
    EuMesh points(2, false);
    for (int in = 0; in < nodes.n_nodes(); ++in) {
        points.add_marker(nodes.coord[in], 0.02);
    }
    auto st_rnk = points.add<int>("rank");
    auto st_idx = points.add<int>("index");
    auto st_min_rank = points.add<int>("min.rank");
    auto st_max_rank = points.add<int>("max.rank");
    for (int in = 0; in < nodes.n_nodes(); ++in) {
        points[in][st_rnk] = nodes.rank[in];
        points[in][st_idx] = nodes.index[in];
        int min_rank = 100000;
        int max_rank = -1000000;
        for (auto inc: nodes.incident.range(in)) {
            int r = nodes.incident.rank[inc];
            min_rank = std::min(min_rank, r);
            max_rank = std::max(max_rank, r);
        }
        points[in][st_min_rank] = min_rank;
        points[in][st_max_rank] = max_rank;
    }
    PvdFile pvd(filename);
    pvd.variables.append("rank", st_rnk);
    pvd.variables.append("index", st_idx);
    pvd.variables.append("min_rank", st_min_rank);
    pvd.variables.append("max_rank", st_max_rank);
    pvd.save(points, 0);
}

int main(int argc, char** argv) {
    mpi::handler handler(argc, argv);
    threads::init(argc, argv);
    threads::info();
    threads::off();

    // Создать сетку
    EuMesh mesh = test9();

    // Сетка с уникальными узлами
    mesh.make_unique_nodes();

    std::cout << "Single process:\n";
    check_mesh(mesh);

    // Файл для записи
    PvdFile pvd("mesh", "output");
    pvd.polyhedral = true;

    auto u = mesh.add<double>("u");
    for (auto cell: mesh) {
        cell[u] = cell.center().norm();
    }

    // Переменные для сохранения
    pvd.variables = {"level", "verts2D"};
    pvd.variables += {"u",  [u](EuCell& cell) -> double { return cell[u]; }};
    pvd.unique_nodes = false;
    pvd.save(mesh, 0.0);

    // Bounding Box для сетки
    Box domain = mesh.bbox();

    // Варианты инициализации ORB декомпозиции
    //ORB orb(domain, "XY", mpi::size());
    //ORB orb(domain, "YX", 13);
    //ORB orb(domain, "YX", 13, 3);
    //orb.use_exact(false);
    RWalk::Ptr orb = RWalk::create(domain, mpi::size());

    // Установить декомпозицию (+ делает redistribute)
    mesh.set_decomposition(orb);

    pvd.save(mesh.locals(), 1.0);

    /*
    PvdFile pvdf("ghosts", "output");
    pvdf.variables = pvd.variables;
    pvdf.save(mesh.ghosts(), 1.0);

    save_markers(mesh.nodes(), "nodes_aft");
    save_markers(mesh.ghost_nodes(), "nodes_aft");
    */

    std::cout << "Distributed:\n";
    check_mesh(mesh);

    return 0;
}
