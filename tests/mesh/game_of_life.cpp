// @brief Клеточный автомат Конвэя, игра "Жизнь".
// Тест демонстрирует использование индексации на структурированных сетках.
#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/pvd_file.h>

using namespace zephyr::io;
using namespace zephyr::geom;
using namespace zephyr::mesh;

using generator::Rectangle;

static Storable<int> u1;
static Storable<int> u2;

void update(EuMesh& mesh) {
   for (int i = 0; i < mesh.nx(); ++i) {
       for (int j = 0; j < mesh.ny(); ++j) {
           auto cell = mesh(i, j);

           // Не очень нравятся такие )( двойные скобки
           int s = mesh(i + 1, j)[u1] +
                   mesh(i - 1, j)[u1] +
                   mesh(i, j + 1)[u1] +
                   mesh(i, j - 1)[u1] +
                   mesh(i + 1, j + 1)[u1] +
                   mesh(i + 1, j - 1)[u1] +
                   mesh(i - 1, j + 1)[u1] +
                   mesh(i - 1, j - 1)[u1];

           if (cell[u1] > 0) {
               cell[u2] = (s == 2 || s == 3) ? 1 : 0;
           }
           else {
               cell[u2] = s == 3 ? 1 : 0;
           }
       }
   }

   for (auto cell: mesh) {
       cell[u1] = cell[u2];
   }
}

void set_random(EuMesh& mesh) {
    for (auto cell: mesh) {
        cell[u1] = rand() % 2;
    }
}

void set_glider(EuMesh& mesh, int i, int j) {
    for (auto cell: mesh) {
        cell[u1] = 0;
    }
    mesh(i + 0, j + 0)[u1] = 1;
    mesh(i + 1, j + 0)[u1] = 1;
    mesh(i + 1, j + 1)[u1] = 1;
    mesh(i + 2, j + 1)[u1] = 1;
    mesh(i + 0, j + 2)[u1] = 1;
}

int main() {
    Rectangle rect(0.0, 1.0, 0.0, 1.0);
    rect.set_nx(50);

    EuMesh mesh(rect);
    u1 = mesh.add<int>("u1");
    u2 = mesh.add<int>("u2");

    PvdFile pvd("life", "output");
    pvd.variables.append("u", u1);

    set_random(mesh);
    //set_glider(mesh, 20, 20);

    int n_steps = 1000;
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