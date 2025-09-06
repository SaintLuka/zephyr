// Тестирование многогранников: построение и сечения

#include <zephyr/geom/primitives/polyhedron.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/pvd_file.h>

using namespace zephyr::geom;
using namespace zephyr::mesh;
using namespace zephyr::io;

void plot_sections(
        const Polyhedron& poly,
        const Vector3d& normal,
        const std::string& figname) {

    PvdFile pvd(figname, "output/" + figname);
    pvd.polyhedral = true;

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

        EuMesh mesh(3, false);
        mesh.push_back(clip);

        pvd.save(mesh, i / (N - 1.0));
    }
}

int main() {
    Vector3d n = {0.1, 0.2, 0.3};
    n.normalize();

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