/// @file Тестирование AMR на задаче с медленно (по времени) и
/// плавно (по координате) изменяющимися уровнями адаптации.

#include <iomanip>

#include <zephyr/mesh/euler/eu_mesh.h>
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
    int wanted;
};

_U_ U;

double get_wanted(AmrStorage::Item& cell) {
    return cell(U).wanted;
}

inline double sqr(double x) {
    return x * x;
}

Vector3d epitrochoid(double t) {
    double R = 0.6;
    double m = 0.2;
    double h = 0.3;
    return {
            R * (m + 1) * std::cos(m * t) - h * std::cos((m + 1) * t),
            R * (m + 1) * std::sin(m * t) - h * std::sin((m + 1) * t),
            0.0
    };
}

// Периодическая функция времени, с периодом = 1
void set_wanted(EuCell& cell, int level, double t) {
    double phi = 10 * M_PI * t;
    Vector3d C = epitrochoid(phi);
    Vector3d T = epitrochoid(phi + 1.0e-6) - epitrochoid(phi - 1.0e-6); // касательная
    Vector3d N = {T.y(), -T.x(), 0.0}; // нормаль
    N.normalize();

    double H = 4.0;

    double dist = std::abs((cell.center() - C).dot(N));

    double xi = std::max(0.0, 1.0 - dist / H);

    int wanted = int(std::floor((level + 0.99) * sqr(sqr(sqr(xi)))));

    cell(U).wanted = wanted;
}

void set_flag(EuCell& cell) {
    int wanted = cell(U).wanted;
    if (cell.level() < wanted) {
        cell.set_flag(1);
    }
    else if (cell.level() == wanted) {
        cell.set_flag(0);
    }
    else {
        cell.set_flag(-1);
    }
}

int main() {
    threads::on();

    PvdFile pvd("mesh", "output");
    pvd.variables = {"index", "level"};
    pvd.variables += { "wanted", get_wanted };

    Rectangle rect(-1.0, 1.0, -1.0, 1.0);
    rect.set_nx(20);

    EuMesh mesh(U, &rect);

    mesh.set_max_level(5);

    int res = mesh.check_base();
    if (res < 0) {
        std::cout << "bad init mesh\n";
        return 0;
    }

    // Начальная адаптация
    std::cout << "Init refinement\n";
    for (int lvl = 0; lvl < mesh.max_level() + 2; ++lvl) {
        std::cout << "  Level " << lvl << "\n";
        mesh.for_each(set_wanted, mesh.max_level(), 0.0);
        mesh.for_each(set_flag);
        mesh.refine();
    }

    Stopwatch elapsed;
    Stopwatch sw_write;
    Stopwatch sw_wanted;
    Stopwatch sw_set_flags;
    Stopwatch sw_refine;

    std::cout << "\nRUN\n";
    elapsed.resume();
    for (int step = 0; step < 1000; ++step) {
        sw_write.resume();
        if (step % 20 == 0) {
            std::cout << "  Step " << std::setw(4) << step << " / 1000\n";
            pvd.save(mesh, step);
        }
        sw_write.stop();

        sw_wanted.resume();
        mesh.for_each(set_wanted, mesh.max_level(), step / 1000.0);
        sw_wanted.stop();

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

    std::cout << "\nElapsed:      " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    std::cout << "  Write:      " << sw_write.extended_time()
              << " ( " << sw_write.milliseconds() << " ms)\n";

    std::cout << "  Set wanted: " << sw_wanted.extended_time()
              << " ( " << sw_wanted.milliseconds() << " ms)\n";

    std::cout << "  Set flags:  " << sw_set_flags.extended_time()
              << " ( " << sw_set_flags.milliseconds() << " ms)\n";

    std::cout << "  Refine:     " << sw_refine.extended_time()
              << " ( " << sw_refine.milliseconds() << " ms)\n";

    return 0;
}