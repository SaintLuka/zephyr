// Тестирование AMR на задаче со случайной адаптацией. В данном тесте флаги
// адаптации выбираются случайным образом с некоторыми вероятностями.

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


// Выставить в ячейке случайный флаг адаптации
void set_flag(EuCell& cell) {
    const double p_coarse = 0.80;  // вероятность огрубления
    const double p_retain = 0.18;  // вероятность сохрнения

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

int main() {
    threads::on();

    Rectangle gen(-1.0, 1.0, -1.0, 1.0);
    gen.set_nx(50);

    EuMesh mesh(gen);
    mesh.set_max_level(4);

    PvdFile pvd("mesh", "output");
    pvd.variables = {"rank", "index", "next", "level", "flag", "faces2D"};

    if (mesh.check_base() < 0) {
        std::cout << "Bad init mesh\n";
        return 0;
    }

    Stopwatch elapsed;
    Stopwatch sw_write;
    Stopwatch sw_set_flags;
    Stopwatch sw_refine;

    std::cout << "RUN\n";
    elapsed.start();
    for (int step = 0; step < 1000; ++step) {
        if (step % 20 == 0) {
            std::cout << "  Step " << std::setw(4) << step << " / 1000\n";
            sw_write.resume();
            pvd.save(mesh, step);
            sw_write.stop();
        }

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

    std::cout << "  Set flags: " << sw_set_flags.extended_time()
              << " ( " << sw_set_flags.milliseconds() << " ms)\n";

    std::cout << "  Refine:    " << sw_refine.extended_time()
              << " ( " << sw_refine.milliseconds() << " ms)\n";

    return 0;
}