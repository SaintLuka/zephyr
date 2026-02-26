/// @brief Тест с полным сохранением/восстановлением сетки

#include <iostream>
#include <iomanip>
#include <random>

#include <zephyr/utils/mpi.h>
#include <zephyr/utils/stopwatch.h>

#include <zephyr/io/pvd_file.h>

#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/geom/generator/rectangle.h>

using namespace zephyr::mesh;
using namespace zephyr::geom;

using zephyr::utils::mpi;
using zephyr::io::PvdFile;
using zephyr::geom::generator::Rectangle;
using zephyr::utils::Stopwatch;

int main(int argc, char** argv) {
    mpi::handler init;

    Rectangle rect(0.0, 1.0, 0.0, 1.0);
    rect.set_nx(20);

    EuMesh mesh(rect);

    auto dens = mesh.add<double>("density");




    return 0;
}