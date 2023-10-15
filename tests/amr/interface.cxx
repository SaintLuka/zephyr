/// @file Тестирование AMR на задачах с интерфейсами, под интерфейсами
/// понимаются узкие области адаптации с шириной в одну ячейку, которые
/// имеют меру размерности меньше, чем мера самой области.

#include <zephyr/mesh/mesh.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/io/pvd_file.h>
#include <zephyr/io/variables.h>
#include <zephyr/utils/stopwatch.h>

using namespace zephyr;
using namespace mesh;

using generator::Rectangle;
using zephyr::io::PvdFile;
using zephyr::utils::Stopwatch;


struct _U_ {
    int idx;
    int bit;
};

_U_ U;

double get_idx(Storage::Item cell) {
    return cell(U).idx;
}

double get_bit(Storage::Item cell) {
    return cell(U).bit;
}

inline double sqr(double x) {
    return x * x;
}

// Периодическая функция времени, с периодом = 1
int calc_idx(ICell& cell, double t) {
    Vector3d c = cell.center();

    double r = c.norm();
    double phi = std::atan2(c.y(), c.x()) - 2 * M_PI * t;
    if (phi < M_PI) {
        phi += M_PI;
    }

    int n_parts = 12;
    double f1 = 0.5 * phi / M_PI + 0.02 * std::sin(12 * r);
    double f2 = 0.5 * phi / M_PI + 0.1 * r * r;
    double sigma = sqr(std::sin(2 * M_PI * t));
    double f3 = sigma * f1 + (1 - sigma) * f2;

    int a1 = (int(std::floor(n_parts * f3)) + n_parts) % n_parts;

    double g1 = std::round(sqr(std::sin(5.0 * r*r - 2*M_PI*t)));
    int a2 = int(g1);
    return a1 + a2;
}

int calc_bit(ICell& cell) {
    int idx_c = cell(U).idx;
    for (auto& face: cell.faces()) {
        if (face.is_boundary()) {
            continue;
        }

        int idx_n = face.neib()(U).idx;
        if (idx_c != idx_n) {
            return 1;
        }
    }
    return 0;
}

void set_index(ICell& cell, double t) {
    cell(U).idx = calc_idx(cell, t);
}

void set_flag(ICell& cell) {
    cell(U).bit = calc_bit(cell);
    cell.set_flag(cell(U).bit > 0 ? 1 : -1);
}

int solution_step(Mesh& mesh, double t = 0.0) {
    for (auto& cell: mesh) {
        if (cell(U).bit > 0) {
            cell.set_flag(1);
        }
        else {
            cell.set_flag(-1);
        }
    }

    mesh.refine();

    for (auto cell: mesh) {
        cell(U).idx = calc_idx(cell, t);
    }
    for (auto cell: mesh) {
        cell(U).bit = calc_bit(cell);
    }

    return mesh.check_refined();
}

int main() {
    threads::off();

    PvdFile pvd("mesh", "output");
    pvd.variables = {"index", "level"};
    pvd.variables += { "idx", get_idx };
    pvd.variables += { "bit", get_bit };

    Rectangle rect(-1.0, 1.0, -1.0, 1.0);
    rect.set_nx(20);

    Mesh mesh(U, &rect);

    mesh.set_max_level(5);

    int res = mesh.check_base();
    if (res < 0) {
        std::cout << "bad init mesh\n";
        return 0;
    }

    // Начальная адаптация
    std::cout << "Начальная адаптация\n";
    for (int lvl = 0; lvl < mesh.max_level() + 2; ++lvl) {
        std::cout << "  Уровень " << lvl << "\n";
        solution_step(mesh);
    }

    Stopwatch elapsed;
    Stopwatch sw_write;
    Stopwatch sw_set_index;
    Stopwatch sw_set_flags;
    Stopwatch sw_refine;

    std::cout << "\nРасчет\n";
    elapsed.resume();
    for (int step = 0; step <= 1000; ++step) {
        sw_write.resume();
        if (step % 10 == 0) {
            std::cout << "  Шаг " << step << " / 1000\n";
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
        //    throw std::runtime_error("Bad mesh");
        //}
    }
    elapsed.stop();

    std::cout << "\nElapsed time:     " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    std::cout << "  Write time:     " << sw_write.extended_time()
              << " ( " << sw_write.milliseconds() << " ms)\n";

    std::cout << "  set index time: " << sw_set_index.extended_time()
              << " ( " << sw_set_index.milliseconds() << " ms)\n";

    std::cout << "  set flags time: " << sw_set_flags.extended_time()
              << " ( " << sw_set_flags.milliseconds() << " ms)\n";

    std::cout << "  Refine time:    " << sw_refine.extended_time()
              << " ( " << sw_refine.milliseconds() << " ms)\n";

    return 0;
}