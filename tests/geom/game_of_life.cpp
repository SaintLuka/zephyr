/// @file Клеточный автомат Конвэя, игра "жизнь".
/// Тест демонстрирует использование индексации для структурированных сеток.

#include <iostream>
#include <iomanip>

#include <zephyr/io/pvd_file.h>
#include <zephyr/mesh/mesh.h>
#include <zephyr/geom/generator/rectangle.h>

using namespace zephyr::io;
using namespace zephyr::geom;
using namespace zephyr::mesh;

using zephyr::mesh::generator::Rectangle;

class _U_ {
public:
    int u1, u2;

    _U_() = default;
};

static const _U_ U{};

void update(EuMesh& mesh) {
   for (int i = 0; i < mesh.nx(); ++i) {
       for (int j = 0; j < mesh.ny(); ++j) {
           auto cell = mesh(i, j);

           // Другой вариант получить данные
           // mesh(i, j)(U).u1
           // но мне не нравятся такие )( двойные скобки
           int s = mesh(i + 1, j).data(U).u1 +
                   mesh(i - 1, j).data(U).u1 +
                   mesh(i, j + 1).data(U).u1 +
                   mesh(i, j - 1).data(U).u1 +
                   mesh(i + 1, j + 1).data(U).u1 +
                   mesh(i + 1, j - 1).data(U).u1 +
                   mesh(i - 1, j + 1).data(U).u1 +
                   mesh(i - 1, j - 1).data(U).u1;

           if (cell(U).u1 > 0) {
               cell(U).u2 = (s == 2 || s == 3) ? 1 : 0;
           }
           else {
               cell(U).u2 = s == 3 ? 1 : 0;
           }
       }
   }

   for (auto cell: mesh) {
       cell(U).u1 = cell(U).u2;
   }
}

int main() {
    PvdFile pvd("life", "out");
    pvd.variables += {"u", [](AmrStorage::Item& cell) -> double { return cell(U).u1; }};

    Rectangle rect(0.0, 1.0, 0.0, 1.0);
    rect.set_nx(50);

    EuMesh mesh(U, &rect);

    for (auto cell: mesh) {
        cell(U).u1 = rand() % 2;
    }

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