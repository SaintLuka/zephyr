// Сохранение массива ячеек и простой эйлеровой сетки.
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/pvd_file.h>
#include <zephyr/utils/mpi.h>

using zephyr::geom::generator::Rectangle;
using zephyr::mesh::EuMesh;
using zephyr::io::Variables;
using zephyr::io::VtuFile;
using zephyr::io::PvdFile;
using zephyr::utils::mpi;

int main(int argc, char** argv) {
    mpi::handler init(argc, argv);

    // Сеточный генератор
    Rectangle gen(0.0, 1.0, 0.0, 1.0, true);
    gen.set_nx(100);

    // Создаем сетку
    EuMesh mesh(gen);
    mesh.set_decomposition("XY");

    // Добавим массив данных к ячейкам
    auto u = mesh.add_cell_data<double>("u");

    // Заполним некоторой функцией
    for (auto cell: mesh) {
        double r = cell.center().norm();
        cell[u] = std::cos(10.0 / (r * r + 0.2));
    }

    // Переменная для записи
    Variables vars = {"faces[6]"};
    vars.add_cell_data("u", u);

    if (mpi::single()) {
        // Базовое сохранение сетки
        VtuFile::save("output/mesh", mesh, vars);

        // Сохранить как набор ячеек
        VtuFile::save("output/cells", mesh.local_cells(), vars);
    }
    else {
        // Сохранение распределенной сетки
        PvdFile pvd("mesh", "output");
        pvd.variables = vars;
        pvd.save(mesh, 0.0);
    }
    return 0;
}