/// @file Тестирование AMR на задачах с интерфейсами, под интерфейсами
/// понимаются узкие области адаптации с шириной в одну ячейку, которые
/// имеют меру размерности меньше, чем мера самой области.

#include <zephyr/mesh/mesh.h>
#include <zephyr/mesh/generator/rectangle.h>
#include <zephyr/io/pvd_file.h>
#include <zephyr/io/variables.h>

using namespace zephyr;
using namespace mesh;

using generator::Rectangle;
using zephyr::io::PvdFile;


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

double get_flag(Storage::Item cell) {
    return cell.flag();
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

int main() {
    PvdFile pvd("mesh", "output");
    pvd.variables += { "idx", get_idx };
    pvd.variables += { "bit", get_bit };
    pvd.variables += { "flag", get_flag };

    Rectangle rect(-1.0, 1.0, -1.0, 1.0);
    rect.set_nx(200);

    Mesh mesh(U, &rect);

    mesh.set_max_level(3);

    int res = mesh.check_base();
    if (res < 0) {
        std::cout << "bad init mesh\n";
        return 0;
    }

    for (int step = 0; step < 2; ++step) {
        std::cout << "Шаг " << step << "\n";
        for (auto cell: mesh.cells()) {
            cell(U).idx = calc_idx(cell, step / 100.0);
        }
        for (auto cell: mesh.cells()) {
            cell(U).bit = calc_bit(cell);
        }
        for (auto& cell: mesh.cells()) {
            if (cell(U).bit > 0) {
                cell.set_flag(1);
            }
            else {
                cell.set_flag(-1);
            }
        }
        pvd.save(mesh, step);
        mesh.refine();
        pvd.save(mesh, step + 0.5);
    }

    return 0;
}