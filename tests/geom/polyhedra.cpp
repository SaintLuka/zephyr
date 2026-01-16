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

    int N = 1000;
    double max_epsilon_cv = 0;
    double avarage_epsilon_cv = 0;
    double max_epsilon_fs = 0;
    double avg_epsilon_fs = 0;
    std::vector<int> checks = std::vector<int>();
    bool flag = false;
    for (int i = 0; i < N; ++i) {
        Vector3d p = (p_max - (p_max - p_min) * i / N) * normal;

        Polyhedron clip = poly.clip(p, normal);

        int check = clip.checkout();
        if (check < 0) {
            std::cout << " p: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
            std::cout << " i: " << i << std::endl;
//            return clip;
            clip = poly.clip(p, normal);
            checks.push_back(check);
//            flag = true;
        }

        double v = poly.volume();
        double v1 = clip.volume();
        double v2 = poly.clip_volume_and_area(p, normal).volume;

        Vector3d p2 = poly.find_section_newton(normal, v1/v);

//        double epsilon_cv = std::abs(v2-v1);
//        avarage_epsilon_cv += epsilon_cv;
//        if (max_epsilon_cv < epsilon_cv)
//            max_epsilon_cv = epsilon_cv;

//        std::cout << " avarage_epsilon_cv: " << avarage_epsilon_cv << std::endl;
//        std::cout << " i: " << i << std::endl;
//        std::cout << std::endl;

        double epsilon_fs = std::abs((p2-p).dot(normal));
        avg_epsilon_fs += epsilon_fs;
        if (max_epsilon_fs < epsilon_fs)
            max_epsilon_fs = epsilon_fs;

//        if (epsilon_fs > 1) {
//            std::cout << " p: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
//            std::cout << " p2: " << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
//            std::cout << " v2/v: " << v1/v << std::endl;
//            std::cout << " i: " << i << std::endl;
//        }

        EuMesh mesh(3, false);
        mesh.push_back(clip);

        pvd.save(mesh, i / (N - 1.0));
        if (flag) {
            std::cout << " p: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
            std::cout << " i: " << i << std::endl;
            break;
        }
    }
//    std::cout << figname << " clip_volume max error: " << max_epsilon_cv << std::endl;
//    std::cout << figname << " clip_volume avg error: " << avarage_epsilon_cv / N << std::endl;
    std::cout << figname << " find_section_newton max error: " << max_epsilon_fs << std::endl;
    std::cout << figname << " find_section_newton avg error: " << avg_epsilon_fs / N << std::endl;
//    std::cout << " checks's size: " << checks.size() << " , checks: ";
    for (size_t i = 0; i < checks.size(); ++i) {
        std::cout << checks[i];
        if (i != checks.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;
//    return poly;
}

int main() {
//    Vector3d n = {0.4, 0.5, 0.3};
//    Vector3d p = {-0.159999, -0.199999, -0.12};
//    Vector3d p = {-0., -3223., -0.4};
//    Vector3d p = {0.16, 0.20, 0.12};
//    n.normalize();
//    auto cube = Polyhedron::Cube();
//    std::cout << Polyhedron::Cube().clip(p, n).checkout()  << std::endl;
//    std::cout << Polyhedron::Cube().center().transpose() << std::endl;

//    Vector3d n = {1, 0, 1};
//    n.normalize();
//    Vector3d p = {0, 0, 0};
//    auto cube = Polyhedron::Cube();
//    std::cout << Polyhedron::Cube().clip(p, n).face_area(Polyhedron::Cube().clip(p, n).n_faces() - 1)  << std::endl;
//    p = {0, 0.001, 0};
//    std::cout << Polyhedron::Cube().clip(p, n).face_area(Polyhedron::Cube().clip(p, n).n_faces() - 1)  << std::endl;

//    Vector3d n = {0., 0., 0.5};
//    n.normalize();
//    std::cout << Polyhedron::Cube().clip_volume(p, n)  << std::endl;
//    std::cout << Polyhedron::Tetrahedron().find_section_newton(n, 0.949853).transpose()  << std::endl;
//    std::cout << Polyhedron::Cube().find_section_newton(n, 0.888431).transpose()  << std::endl;

    Vector3d n = {0.4, 0.5, 0.3};
    n.normalize();
    plot_sections(Polyhedron::Cube(), n, "cube");
//    plot_sections(Polyhedron::Cuboid(1.0, 1.2, 0.8), n, "cuboid");
//    plot_sections(Polyhedron::Pyramid(), n, "pyramid");
//    plot_sections(Polyhedron::Wedge(), n, "wedge");
//    plot_sections(Polyhedron::Tetrahedron(), n, "tetrahedron");
//    plot_sections(Polyhedron::Octahedron(), n, "octahedron");
//    plot_sections(Polyhedron::Dodecahedron(), n, "dodecahedron");
//    plot_sections(Polyhedron::Icosahedron(), n, "icosahedron");
//    plot_sections(Polyhedron::TruncatedCube(), n, "truncated_cube");
//
//    n = {-1, 0, 1};
//    n.normalize();
//
//    plot_sections(Polyhedron::Cube(), n, "cube");
//    plot_sections(Polyhedron::Cuboid(1.0, 1.2, 0.8), n, "cuboid");
//    plot_sections(Polyhedron::Pyramid(), n, "pyramid");
//    plot_sections(Polyhedron::Wedge(), n, "wedge");
//    plot_sections(Polyhedron::Tetrahedron(), n, "tetrahedron");
//    plot_sections(Polyhedron::Octahedron(), n, "octahedron");
//    plot_sections(Polyhedron::Dodecahedron(), n, "dodecahedron");
//    plot_sections(Polyhedron::Icosahedron(), n, "icosahedron");
//    plot_sections(Polyhedron::TruncatedCube(), n, "truncated_cube");
//
//    n = {0.5, 0.5, 0.5};
//    n.normalize();
//
//    plot_sections(Polyhedron::Cube(), n, "cube");
//    plot_sections(Polyhedron::Cuboid(1.0, 1.2, 0.8), n, "cuboid");
//    plot_sections(Polyhedron::Pyramid(), n, "pyramid");
//    plot_sections(Polyhedron::Wedge(), n, "wedge");
//    plot_sections(Polyhedron::Tetrahedron(), n, "tetrahedron");
//    plot_sections(Polyhedron::Octahedron(), n, "octahedron");
//    plot_sections(Polyhedron::Dodecahedron(), n, "dodecahedron");
//    plot_sections(Polyhedron::Icosahedron(), n, "icosahedron");
//    plot_sections(Polyhedron::TruncatedCube(), n, "truncated_cube");
//
//    n = {0.5, 0, 0};
//    n.normalize();
//
//    plot_sections(Polyhedron::Cube(), n, "cube");
//    plot_sections(Polyhedron::Cuboid(1.0, 1.2, 0.8), n, "cuboid");
//    plot_sections(Polyhedron::Pyramid(), n, "pyramid");
//    plot_sections(Polyhedron::Wedge(), n, "wedge");
//    plot_sections(Polyhedron::Tetrahedron(), n, "tetrahedron");
//    plot_sections(Polyhedron::Octahedron(), n, "octahedron");
//    plot_sections(Polyhedron::Dodecahedron(), n, "dodecahedron");
//    plot_sections(Polyhedron::Icosahedron(), n, "icosahedron");
//    plot_sections(Polyhedron::TruncatedCube(), n, "truncated_cube");
//
//    n = {0.3, 12, 53};
//    n.normalize();
//
//    plot_sections(Polyhedron::Cube(), n, "cube");
//    plot_sections(Polyhedron::Cuboid(1.0, 1.2, 0.8), n, "cuboid");
//    plot_sections(Polyhedron::Pyramid(), n, "pyramid");
//    plot_sections(Polyhedron::Wedge(), n, "wedge");
//    plot_sections(Polyhedron::Tetrahedron(), n, "tetrahedron");
//    plot_sections(Polyhedron::Octahedron(), n, "octahedron");
//    plot_sections(Polyhedron::Dodecahedron(), n, "dodecahedron");
//    plot_sections(Polyhedron::Icosahedron(), n, "icosahedron");
//    plot_sections(Polyhedron::TruncatedCube(), n, "truncated_cube");

    return 0;
}