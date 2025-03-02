/// @file Тестирование AMR на задачах с движением больших областей с высоким
/// уровнем адаптации. Под большими областями понимаются области, которые имеют
/// ту же меру размерности, что и сама область.

#include <iomanip>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/polygon.h>
#include <zephyr/io/pvd_file.h>
#include <zephyr/io/variables.h>
#include <zephyr/utils/stopwatch.h>

using namespace zephyr;
using namespace mesh;

using generator::Rectangle;
using zephyr::io::PvdFile;
using zephyr::utils::Stopwatch;


struct _U_ {
    int bit;
};

_U_ U;

double get_bit(AmrStorage::Item& cell) {
    return cell(U).bit;
}

inline double sqr(double x) {
    return x * x;
}

/// @brief Пятиконечная звезда
struct Star {
    double R;
    Polygon poly;

    /// @param R Радиус описаной окружности
    explicit Star(double R) : R(R), poly(10) {
        update(0.0);
    }

    /// @brief Точка внутри полигона?
    bool inside(const Vector3d& v) const {
        return poly.inside(v);
    }

    void update(double t) {
        double xi = 2.0 * M_PI * t;
        double phi = 0.00123 + xi;

        double r2 = R;
        double r1 = 0.4 * R;
        double r = r1 + (r2 - r1) * sqr(std::sin(0.5 * xi));

        double a = 0.9;
        double xc = a * std::sqrt(2.0) * std::cos(xi) / (1.0 + sqr(std::sin(xi)));
        double yc = xc * std::sin(xi);

        for (int i = 0; i < 10; i += 2) {
            Vector3d v1 = {xc + R * std::cos(M_PI * i / 5.0 + phi),
                           yc + R * std::sin(M_PI * i / 5.0 + phi), 0.0};

            Vector3d v2 = {xc + r * std::cos(M_PI * (i + 1) / 5.0 + phi),
                           yc + r * std::sin(M_PI * (i + 1) / 5.0 + phi), 0.0};

            poly.set(i, v1);
            poly.set(i + 1, v2);
        }
    }
};

void set_index(EuCell& cell, Star& star) {
    cell(U).bit = star.inside(cell.center());
}

void set_flag(EuCell& cell) {
    cell.set_flag(cell(U).bit > 0 ? 1 : -1);
}

int main() {
    threads::on();

    PvdFile pvd("mesh", "output");
    pvd.variables = {"index", "level"};
    pvd.variables += { "bit", get_bit };

    Rectangle rect(-2.0, 2.0, -1.0, 1.0);
    rect.set_nx(50);

    EuMesh mesh(rect, U);

    mesh.set_max_level(5);

    int res = mesh.check_base();
    if (res < 0) {
        std::cout << "bad init mesh\n";
        return 0;
    }

    // Фигура для адаптации
    Star star(0.5);

    // Начальная адаптация
    std::cout << "Init refinement\n";
    for (int lvl = 0; lvl < mesh.max_level() + 4; ++lvl) {
        std::cout << "  Level " << lvl << "\n";
        mesh.for_each(set_index, star);
        mesh.for_each(set_flag);
        mesh.refine();
    }

    Stopwatch elapsed;
    Stopwatch sw_write;
    Stopwatch sw_star;
    Stopwatch sw_set_index;
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

        sw_star.resume();
        star.update(step / 1000.0);
        sw_star.stop();

        sw_set_index.resume();
        mesh.for_each(set_index, star);
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

    std::cout << "\nElapsed:       " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    std::cout << "  Write:       " << sw_write.extended_time()
              << " ( " << sw_write.milliseconds() << " ms)\n";

    std::cout << "  Update star: " << sw_star.extended_time()
              << " ( " << sw_star.milliseconds() << " ms)\n";

    std::cout << "  Set index:   " << sw_set_index.extended_time()
              << " ( " << sw_set_index.milliseconds() << " ms)\n";

    std::cout << "  Set flags:   " << sw_set_flags.extended_time()
              << " ( " << sw_set_flags.milliseconds() << " ms)\n";

    std::cout << "  Refine:      " << sw_refine.extended_time()
              << " ( " << sw_refine.milliseconds() << " ms)\n";

    threads::off();
    return 0;
}