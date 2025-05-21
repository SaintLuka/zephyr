/// @file Тестирование AMR на задаче с медленно (по времени) и
/// плавно (по координате) изменяющимися уровнями адаптации.

#include <iomanip>

#include <zephyr/mesh/euler/soa_mesh.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/io/pvd_file.h>
#include <zephyr/utils/stopwatch.h>

using namespace zephyr;
using namespace mesh;

using generator::Rectangle;
using zephyr::io::PvdFile;
using zephyr::utils::Stopwatch;

// Поле для хранения на сетке
static Storable<int> wanted;

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
void set_wanted(QCell& cell, int level, double t) {
    double phi = 10 * M_PI * t;
    Vector3d C = epitrochoid(phi);
    Vector3d T = epitrochoid(phi + 1.0e-6) - epitrochoid(phi - 1.0e-6); // касательная
    Vector3d N = {T.y(), -T.x(), 0.0}; // нормаль
    N.normalize();

    double H = 4.0;

    double dist = std::abs((cell.center() - C).dot(N));

    double xi = std::max(0.0, 1.0 - dist / H);

    cell(wanted) = int(std::floor((level + 0.99) * std::pow(xi, 8)));
}

void set_flag(QCell& cell) {
    if (cell.level() < cell(wanted)) {
        cell.set_flag(1);
    }
    else if (cell.level() == cell(wanted)) {
        cell.set_flag(0);
    }
    else {
        cell.set_flag(-1);
    }
}

int main() {
    threads::on();

    PvdFile pvd("mesh", "output");
    pvd.variables = {"rank", "index", "next", "level", "flag", "faces2D"};
    pvd.variables += {"wanted", [](QCell& cell) -> double { return cell(wanted); }};

    Rectangle rect(-1.0, 1.0, -1.0, 1.0);
    rect.set_nx(20);

    SoaMesh mesh(rect);
    wanted = mesh.add_data<int>("wanted");

    mesh.set_max_level(5);

    if (mesh.check_base() < 0) {
        std::cout << "Bad init mesh\n";
        return 0;
    }

    // Начальная адаптация
    std::cout << "Init refinement\n";
    for (int lvl = 0; lvl < mesh.max_level() + 2; ++lvl) {
        std::cout << "  Level " << lvl << "\n";
        mesh.for_each(set_wanted, mesh.max_level(), 0.0);
        mesh.for_each(set_flag);
        mesh.refine();

        //if (mesh.check_refined() < 0) {
        //    throw std::runtime_error("Bad init refinement");
        //}
    }

    Stopwatch elapsed;
    Stopwatch sw_write;
    Stopwatch sw_wanted;
    Stopwatch sw_set_flags;
    Stopwatch sw_refine;

    std::cout << "RUN\n";
    elapsed.start();
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
        //    throw std::runtime_error("Bad refined mesh");
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