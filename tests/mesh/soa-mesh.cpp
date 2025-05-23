#include <zephyr/utils/threads.h>
#include <zephyr/mesh/euler/soa_mesh.h>
#include <zephyr/io/vtu_file.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/geom/generator/cuboid.h>

using zephyr::utils::threads;
using namespace zephyr::io;
using zephyr::mesh::SoaMesh;
using zephyr::mesh::SoaStorage;
using zephyr::mesh::AmrCells;
using zephyr::geom::generator::Rectangle;
using zephyr::geom::generator::Cuboid;

using zephyr::geom::Boundary;


int main() {
    threads::on();

    //Cuboid gen(0.0, 1.0, 0.0, 0.6, 0.0, 0.9);
    Rectangle gen(-0.5, 1.0, -0.3, 0.7);
    gen.set_sizes(100, 100);
    gen.set_boundaries({.left=Boundary::ZOE,.right=Boundary::ZOE,
                        .bottom=Boundary::ZOE,.top=Boundary::WALL});

    SoaMesh mesh(gen);
    mesh.add<double>("value");

    if (mesh.check_base() < 0) {
        std::cout << "Bad mesh.\n";
    }
    VtuFile vtu("mesh.vtu");
    vtu.variables = {"faces2D", "rank", "level", "index", "flag", "b_idx", "z_idx"};
    vtu.save(mesh);

    return 0;
}