/// @file Тестирование AMR на задаче с медленно (по времени) и
/// плавно (по координате) изменяющимися уровнями адаптации.

#include <zephyr/mesh/mesh.h>
#include <zephyr/mesh/generator/rectangle.h>
#include <zephyr/io/pvd_file.h>
#include <zephyr/io/variables.h>

using namespace zephyr;
using namespace mesh;

using generator::Rectangle;
using zephyr::io::PvdFile;


struct _U_ {
    int wanted;
};

_U_ U;

double get_wanted(Storage::Item cell) {
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
int calc_wanted(ICell& cell, int level, double t) {
    double phi = 10 * M_PI * t;
    Vector3d C = epitrochoid(phi);
    Vector3d T = epitrochoid(phi + 1.0e-6) - epitrochoid(phi - 1.0e-6); // касательная
    Vector3d N = {T.y(), -T.x(), 0.0}; // нормаль
    N.normalize();

    double H = 4.0;

    double dist = std::abs((cell.center() - C).dot(N));

    double xi = std::max(0.0, 1.0 - dist / H);

    int wanted = int(std::floor((level + 0.99) * sqr(sqr(sqr(xi)))));

    return std::max(0, std::min(wanted, level));
}

int solution_step(Mesh& mesh, double t = 0.0) {
    for (auto& cell: mesh.cells()) {
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

    mesh.refine();

    for (auto cell: mesh.cells()) {
        cell(U).wanted = calc_wanted(cell, mesh.max_level(), t);
    }

    return mesh.check_refined();
}

int main() {
    PvdFile pvd("mesh", "output");
    pvd.variables = {"index", "level"};
    pvd.variables += { "wanted", get_wanted };

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
    std::cout << "Начальная адапция\n";
    for (int lvl = 0; lvl < mesh.max_level() + 2; ++lvl) {
        std::cout << "  Уровень " << lvl << "\n";
        solution_step(mesh);
    }

    std::cout << "\nРасчет\n";
    for (int step = 0; step < 1000; ++step) {
        if (step % 10 == 0) {
            std::cout << "  Шаг " << step << " / 1000\n";
            pvd.save(mesh, step);
        }

        solution_step(mesh, step / 1000.0);
    }

    return 0;
}