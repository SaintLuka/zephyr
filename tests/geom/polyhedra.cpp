/// @brief Тестирование многогранников, сечения, объемы и прочее

#include <zephyr/geom/primitives/polyhedron.h>
#include <zephyr/io/pvd_file.h>

using namespace zephyr::geom;
using namespace zephyr::io;

void plot_sections(
        const Polyhedron& poly,
        const Vector3d& normal,
        const std::string& figname) {

    PvdFile pvd(figname, "output/" + figname);
    pvd.polyhedral = true;

    AmrStorage cell(1);

    double p_min = +1.0e100;
    double p_max = -1.0e100;
    for (int i = 0; i < poly.n_verts(); ++i) {
        p_min = std::min(p_min, poly.vertex(i).dot(normal));
        p_max = std::max(p_max, poly.vertex(i).dot(normal));
    }
    p_min -= 0.001 * (p_max - p_min);
    p_max += 0.001 * (p_max - p_min);

    int N = 101;
    for (int i = 0; i < N; ++i) {
        Vector3d p = (p_max - (p_max - p_min) * i / N) * normal;

        Polyhedron clip = poly.clip(p, normal);

        cell[0] = zephyr::mesh::AmrCell(clip);
        pvd.save(cell, i);
    }
}

int main() {
    //Vector3d n = {0., 0., 0.5};
    //Vector3d point = {0.25, 0.25, 0.5};
    Vector3d n = {0.4, 0.5, 0.1};
    Vector3d point = {0.4, 0.4, 0.2};
    n.normalize();

    Polyhedron cube = Polyhedron::Cube();
    std::cout << cube.clip_volume(point, n) << std::endl;
    std::cout << cube.clip(point, n).volume() << std::endl;

    plot_sections(Polyhedron::Cube(), n, "cube");
    plot_sections(Polyhedron::Cuboid(1.0, 1.2, 0.8), n, "cuboid");
    plot_sections(Polyhedron::Pyramid(), n, "pyramid");
    plot_sections(Polyhedron::Wedge(), n, "wedge");
    plot_sections(Polyhedron::Tetrahedron(), n, "tetrahedron");
    plot_sections(Polyhedron::Octahedron(), n, "octahedron");
    plot_sections(Polyhedron::Dodecahedron(), n, "dodecahedron");
    plot_sections(Polyhedron::Icosahedron(), n, "icosahedron");
    plot_sections(Polyhedron::TruncatedCube(), n, "truncated_cube");

    return 0;
}