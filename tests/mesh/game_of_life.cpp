// @brief Клеточный автомат Конвэя, игра "Жизнь".
// Тест демонстрирует использование индексации на структурированных сетках.
#include <iostream>
#include <iomanip>
#include <set>

#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/pvd_file.h>

using namespace zephyr::io;
using namespace zephyr::geom;
using namespace zephyr::mesh;
using namespace zephyr::utils;

using generator::Rectangle;

static Storable<int> u1;
static Storable<int> u2;

void update(EuMesh& mesh) {
    mesh.sync_cells(u1);
    mesh.for_each([](EuCell& cell) {
        int s = cell.neib(+1,  0)[u1] +
                cell.neib(-1,  0)[u1] +
                cell.neib( 0, +1)[u1] +
                cell.neib( 0, -1)[u1] +
                cell.neib(+1, +1)[u1] +
                cell.neib(+1, -1)[u1] +
                cell.neib(-1, +1)[u1] +
                cell.neib(-1, -1)[u1];

        if (cell[u1] > 0) {
            cell[u2] = (s == 2 || s == 3) ? 1 : 0;
        }
        else {
            cell[u2] = s == 3 ? 1 : 0;
        }
    });
    mesh.swap(u1, u2);
}

void set_zero(EuMesh& mesh) {
    for (auto cell: mesh) {
        cell[u1] = 0;
    }
}

void set_random(EuMesh& mesh) {
    for (auto cell: mesh) {
        cell[u1] = rand() % 2;
    }
}

using field_t = std::set<std::tuple<int, int>>;

void add_glider(field_t& field, int i, int j, bool inv = false) {
    int sgn = inv ? -1 : +1;
    field.insert({i-sgn, j});
    field.insert({i,     j});
    field.insert({i,     j+1});
    field.insert({i+sgn, j+1});
    field.insert({i-sgn, j+2});
}

void add_spaceship(field_t& field, int i, int j, bool inv = false) {
    int sgn = inv ? -1 : +1;
    field.insert({i-2*sgn, j + 1});
    field.insert({i-2*sgn, j + 3});
    field.insert({i-1*sgn, j + 0});
    field.insert({i+0*sgn, j + 0});
    field.insert({i+1*sgn, j + 0});
    field.insert({i+2*sgn, j + 0});
    field.insert({i+2*sgn, j + 1});
    field.insert({i+2*sgn, j + 2});
    field.insert({i+1*sgn, j + 3});
}

void set_flottila(EuMesh& mesh) {
    field_t field;
    for (int k = 1; k < 7; ++k) {
        int i1 = 35 + 2 * k;
        int i2 = 65 - 2 * k;
        int j = 15 * k - 1;        
        add_glider(field, i1, j, true);
        add_glider(field, i2, j, false);
    }
    for (int k = 1; k < 7; ++k) {
        int i1 = 3 + 2 * k;
        int i2 = 97 - 2 * k;
        int j = 15 * k - 5;
        add_spaceship(field, i1, j, false);
        add_spaceship(field, i2, j, true);
    }

    for (auto cell: mesh) {
        int i = std::floor(cell.x());
        int j = std::floor(cell.y());
        cell[u1] = field.contains({i, j});
    }
}

int main(int argc, char** argv) {
    mpi::handler handler(argc, argv);
    threads::init(argc, argv);
    threads::info();
    threads::off();

    Rectangle rect(0.0, 100.0, 0.0, 100.0);
    rect.set_nx(100);

    EuMesh mesh(rect, true);
    mesh.set_decomposition("XY");

    u1 = mesh.add<int>("u1");
    u2 = mesh.add<int>("u2");

    PvdFile pvd("life", "output");
    pvd.variables.add_cell_data("u", u1);

    //set_random(mesh);
    set_flottila(mesh);

    int n_steps = 300;
    for (int step = 0; step < n_steps; ++step) {
        if (step % 20 == 0) {
            std::cout << "  Step " << std::setw(4) << step << " / " << n_steps << "\n";
        }
        pvd.save(mesh, step);
        update(mesh);
    }
    pvd.save(mesh, n_steps);

    return 0;
}