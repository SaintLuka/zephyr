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
    double max_epsilon_fs = 0;
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
        double v2 = poly.clip_volume(p, normal);

        Vector3d p2 = poly.find_section(normal, v2/v);

        double epsilon_cv = std::abs(v2-v1);
        if (max_epsilon_cv < epsilon_cv)
            max_epsilon_cv = epsilon_cv;

        double epsilon_fs = std::abs((p2-p).dot(normal)); //???
        if (max_epsilon_fs < epsilon_fs)
            max_epsilon_fs = epsilon_fs;


        EuMesh mesh(3, false);
        mesh.push_back(clip);

        pvd.save(mesh, i / (N - 1.0));
        if (flag) {
            std::cout << " p: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
            std::cout << " i: " << i << std::endl;
            break;
        }
    }
//    std::cout << figname << " clip_volume error: " << max_epsilon_cv << std::endl;
//    std::cout << figname << " find_section error: " << max_epsilon_fs << std::endl;
    std::cout << " checks's size: " << checks.size() << " , checks: ";
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
    Vector3d n = {0., 0.3, 0.3}; //{0.4, 0.5, 0.3};
//    Vector3d p = {-0.159999, -0.199999, -0.12};
//    Vector3d p = {0.16, 0.20, 0.12};
    n.normalize();
//    auto cube = Polyhedron::Cube();
//    std::cout << Polyhedron::Cube().clip(p, n).checkout()  << std::endl;

//    auto wrong_clip = plot_sections(Polyhedron::Cube(), n, "cube");
    plot_sections(Polyhedron::Cube(), n, "cube");
    plot_sections(Polyhedron::Cuboid(1.0, 1.2, 0.8), n, "cuboid");
    plot_sections(Polyhedron::Pyramid(), n, "pyramid");
    plot_sections(Polyhedron::Wedge(), n, "wedge");
    plot_sections(Polyhedron::Tetrahedron(), n, "tetrahedron");
    plot_sections(Polyhedron::Octahedron(), n, "octahedron");
    plot_sections(Polyhedron::Dodecahedron(), n, "dodecahedron");
    plot_sections(Polyhedron::Icosahedron(), n, "icosahedron");
    plot_sections(Polyhedron::TruncatedCube(), n, "truncated_cube");

    n = {0.4, 0.5, 0.3};
//    Vector3d p = {-0.159999, -0.199999, -0.12};
//    Vector3d p = {0.16, 0.20, 0.12};
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

//    std::cout << "Number of vertices: " << wrong_clip.n_verts() << std::endl;
//    std::cout << "Number of faces: " << wrong_clip.n_faces() << std::endl;
//    std::cout << "Volume: " << wrong_clip.volume() << std::endl;
//
//    std::cout << "Center: " << wrong_clip.center().transpose() << std::endl;
//
//    int check_result = wrong_clip.checkout();
//    std::cout << "Checkout result: " << check_result << std::endl;
//
//    std::cout << "\n--- Vertices ---" << std::endl;
//    for (int i = 0; i < wrong_clip.n_verts(); ++i) {
//        std::cout << "Vertex " << i << ": " << wrong_clip.vertex(i).transpose() << std::endl;
//    }
//
//    std::cout << "\n--- Faces ---" << std::endl;
//    for (int i = 0; i < wrong_clip.n_faces(); ++i) {
//        const std::vector<int>& face_indices = wrong_clip.face_indices(i);
//        std::cout << "Face " << i << " (size " << face_indices.size() << "): ";
//        for (int idx : face_indices) {
//            std::cout << idx << " ";
//        }
//        std::cout << std::endl;
//        std::cout << "  Center: " << wrong_clip.face_center(i).transpose() << std::endl;
//        std::cout << "  Area: " << wrong_clip.face_area(i) << std::endl;
//        std::cout << "  Normal: " << wrong_clip.face_normal(i).transpose() << std::endl;
//        std::cout << "  Area vector: " << wrong_clip.face_area_n(i).transpose() << std::endl;
//    }

//    Vector3d norm = {1.0, 1.0, 1.0};
//    Vector3d p = {0.00000000001, 0., 0.};
//
//    Polyhedron cube = Polyhedron::Cuboid(2.0, 2.0, 2.0);
//    Polyhedron cliped_cube = cube.clip(p, norm);
//    std::cout << cube.checkout() << std::endl;
//    std::cout << cube.volume() << std::endl;


    /*Polyhedron cube = Polyhedron::Cuboid(2.0, 2.0, 2.0);
    Vector3d section_point_c = cube.find_section(norm, 0.5);
    std::cout << "cube.clip_volume(find_section):" << cube.clip_volume(section_point_c, norm) << std::endl;
    std::cout << "cube.clip_volume(p):" << cube.clip_volume(p, norm) << std::endl;
    std::cout << "cube.clip(find_section).volume():" << cube.clip(section_point_c, norm).volume() << std::endl;
    std::cout << "cube.clip(p).volume():" << cube.clip(p, norm).volume() << std::endl;

    Polyhedron pyramid = Polyhedron::Pyramid();
    Vector3d section_point_p = pyramid.find_section(norm, 0.5);
    std::cout << "pyramid.clip_volume(find_section):" << pyramid.clip_volume(section_point_p, norm) << std::endl;
    std::cout << "pyramid.clip_volume(p):" << pyramid.clip_volume(p, norm) << std::endl;
    //std::cout << "pyramid.clip(find_section).volume():" << pyramid.clip(section_point_p, norm).volume() << std::endl;
    //std::cout << "pyramid.clip(p).volume():" << pyramid.clip(p, norm).volume() << std::endl;

    Polyhedron icosahedron = Polyhedron::Icosahedron();
    Vector3d section_point_i = icosahedron.find_section(norm, 0.5);
    std::cout << "icosahedron.clip_volume(find_section):" << icosahedron.clip_volume(section_point_i, norm) << std::endl;
    std::cout << "icosahedron.clip_volume(p):" << icosahedron.clip_volume(p, norm) << std::endl;
    //std::cout << "icosahedron.clip(find_section).volume():" << icosahedron.clip(section_point_i, norm).volume() << std::endl;
    //std::cout << "icosahedron.clip(p).volume():" << icosahedron.clip(p, norm).volume() << std::endl;
    */
    return 0;
}