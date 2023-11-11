/// @file Тестирование AMR на задаче со случайной адаптацией.
/// В данном тесте флаги адаптации выбираются случайным образом
/// с некоторыми вероятностями.

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/io/pvd_file.h>
#include <zephyr/io/variables.h>

using namespace zephyr;
using namespace mesh;

using generator::Rectangle;
using generator::Cuboid;
using zephyr::io::PvdFile;


struct _U_ {
    int val;
};

_U_ U;

double get_flag(AmrStorage::Item& cell) {
    return cell.flag;
}

double get_next(AmrStorage::Item& cell) {
    return cell.next;
}

double get_some(AmrStorage::Item& cell) {
    return cell.index;
}

inline double sqr(double x) {
    return x * x;
}

int solution_step(EuMesh& mesh) {
    const double p_coarse = 0.80;
    const double p_retain = 0.18;
    const double p_refine = 0.02;

    for (auto& cell: mesh) {
        double p = rand() / double(RAND_MAX);

        if (p < p_coarse) {
            cell.set_flag(-1);
        }
        else if (p < p_coarse + p_retain) {
            cell.set_flag(0);
        }
        else {
            cell.set_flag(+1);
        }
    }

    mesh.refine();

    return mesh.check_refined();
}

int main() {
    threads::on();

    PvdFile pvd("mesh", "output");
    pvd.variables = {"index", "level"};

    Rectangle rect(-1.0, 1.0, -1.0, 1.0);
    rect.set_nx(20);

    Cuboid cube(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    cube.set_nx(10);

    EuMesh mesh(U, &cube);

    mesh.set_max_level(3);

    int res = mesh.check_base();
    if (res < 0) {
        std::cout << "bad init mesh\n";
        return 0;
    }

    for (int step = 0; step < 200; ++step) {
        std::cout << "  Шаг " << step << " / 200\n";
        pvd.save(mesh, step);

        solution_step(mesh);
    }

    return 0;
}