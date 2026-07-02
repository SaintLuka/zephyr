// Проверяет average_flux, сравнивает с сечениями
#include <zephyr/utils/tests.h>
#include <zephyr/geom/geom.h>
#include <zephyr/geom/sections.h>
#include <zephyr/geom/primitives/polyhedron.h>

using namespace zephyr::geom;

TEST_CASE("average_flux 2D") {

}

TEST_CASE("average_flux 3D to 2D reduce") {

}

TEST_CASE("average_flux 3D") {

    Polyhedron cube = Polyhedron::Cube();
    //cube.move({0.1, 0.23, 0.77});

    std::cout << cube.volume() << "\n";

    Vector3d cell_c = cube.center();

    double alpha = 0.1;
    Vector3d ni = {0.2, 0.4, 0.3};
    ni.normalize();

    Vector3d p = cell_c + cube_find_section(alpha, ni) * ni;


    double CFL = 1.0e-9;


    for (int i = 0; i < 6; ++i) {
        Vector3d fc = cube.face_center(i);
        Vector3d face_n = cube.face_normal(i);

        Polyhedron clip = cube.clip(fc + 2 * CFL * (cell_c - fc), -face_n);

        std::cout << "CFL: " << CFL << " " << clip.volume() << "\n";

        std::cout << clip.clip(p, ni).volume() / clip.volume() << " " << average_flux(alpha, ni, face_n, CFL) << "\n";
    }


}

int main() {
    return TestRegistry::instance().run_all();
}