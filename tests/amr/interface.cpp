// Тестирование AMR на задачах с интерфейсами, под интерфейсами понимаются
// узкие области адаптации с шириной в одну ячейку, которые имеют меру
// размерности меньше, чем мера самой области.

#include <iomanip>

#include <zephyr/utils/stopwatch.h>
#include <zephyr/utils/threads.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/pvd_file.h>

using namespace zephyr::mesh;
using namespace zephyr::geom;
using namespace zephyr::utils;

using generator::Rectangle;
using zephyr::io::PvdFile;


// Поля для хранения на сетке
static Storable<int> idx;
static Storable<int> bit;

// Периодическая функция времени, с периодом = 1
int calc_idx(EuCell& cell, double t) {
    Vector3d c = cell.center();

    double r = c.norm();
    double phi = std::atan2(c.y(), c.x()) - 2 * M_PI * t;
    if (phi < M_PI) {
        phi += M_PI;
    }

    int n_parts = 12;
    double f1 = 0.5 * phi / M_PI + 0.02 * std::sin(12 * r);
    double f2 = 0.5 * phi / M_PI + 0.1 * r * r;
    double sigma = std::pow(std::sin(2 * M_PI * t), 2);
    double f3 = sigma * f1 + (1 - sigma) * f2;

    int a1 = (int(std::floor(n_parts * f3)) + n_parts) % n_parts;

    double g1 = std::round(std::pow(std::sin(5.0 * r*r - 2*M_PI*t), 2));
    int a2 = int(g1);
    return a1 + a2;
}

// Отметить ячейки на границах областей
int calc_bit(EuCell& cell) {
    int idx_c = cell(idx);
    for (auto& face: cell.faces()) {
        if (face.is_boundary()) {
            continue;
        }

        int idx_n = face.neib()(idx);
        if (idx_c != idx_n) {
            return 1;
        }
    }
    return 0;
}

// Выставить функцию-индикатор подобласти
void set_index(EuCell& cell, double t) {
    cell(idx) = calc_idx(cell, t);
}

void set_flag(EuCell& cell) {
    cell(bit) = calc_bit(cell);
    cell.set_flag(cell(bit) > 0 ? 1 : -1);
}

int main() {
    threads::on();

    Rectangle rect(-1.0, 1.0, -1.0, 1.0);
    rect.set_nx(30);

    EuMesh mesh(rect);
    idx = mesh.add<int>("idx");
    bit = mesh.add<int>("bit");

    mesh.set_max_level(4);
    mesh.set_distributor("simple");

    PvdFile pvd("mesh", "output");
    pvd.variables = {"rank", "index", "next", "level", "flag", "faces2D"};
    pvd.variables += {"idx", [](EuCell& cell) -> double { return cell(idx); }};
    pvd.variables += {"bit", [](EuCell& cell) -> double { return cell(bit); }};

    if (mesh.check_base() < 0) {
        std::cout << "Bad init mesh\n";
        return 0;
    }

    // Начальная адаптация
    std::cout << "Init refinement\n";
    for (int lvl = 0; lvl < mesh.max_level() + 2; ++lvl) {
        std::cout << "  Level " << lvl << "\n";

        mesh.for_each(set_index, 0.0);
        mesh.for_each(set_flag);
        mesh.refine();

        if (mesh.check_refined() < 0) {
            std::cout << "Bad init refinement\n";
            return 0;
        }
    }

    Stopwatch elapsed;
    Stopwatch sw_write;
    Stopwatch sw_set_index;
    Stopwatch sw_set_flags;
    Stopwatch sw_refine;

    std::cout << "RUN\n";
    elapsed.start();
    for (int step = 0; step <= 1000; ++step) {
        sw_write.resume();
        if (step % 20 == 0) {
            std::cout << "  Step " << std::setw(4) << step << " / 1000\n";
            pvd.save(mesh, step);
        }
        sw_write.stop();

        sw_set_index.resume();
        mesh.for_each(set_index, step / 1000.0);
        sw_set_index.stop();

        sw_set_flags.resume();
        mesh.for_each(set_flag);
        sw_set_flags.stop();

        sw_refine.resume();
        mesh.refine();
        sw_refine.stop();

        //if (mesh.check_refined() < 0) {
        //    throw std::runtime_error("Bad refined mesh");
        //}
    }
    elapsed.stop();

    std::cout << "\nElapsed:     " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    std::cout << "  Write:     " << sw_write.extended_time()
              << " ( " << sw_write.milliseconds() << " ms)\n";

    std::cout << "  Set index: " << sw_set_index.extended_time()
              << " ( " << sw_set_index.milliseconds() << " ms)\n";

    std::cout << "  Set flags: " << sw_set_flags.extended_time()
              << " ( " << sw_set_flags.milliseconds() << " ms)\n";

    std::cout << "  Refine:    " << sw_refine.extended_time()
              << " ( " << sw_refine.milliseconds() << " ms)\n";

    return 0;
}