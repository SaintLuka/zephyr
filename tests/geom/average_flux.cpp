// Проверяет average_flux, сравнивает с сечениями
#include <zephyr/utils/tests.h>
#include <zephyr/geom/geom.h>
#include <zephyr/geom/sections.h>
#include <zephyr/geom/primitives/polyhedron.h>

using namespace zephyr::geom;

TEST_CASE("average_flux 2D") {
    std::array vs = {
        Vector3d{-1, -1, 0},
        Vector3d{+1, -1, 0},
        Vector3d{+1, +1, 0},
        Vector3d{-1, +1, 0}
    };
    for (auto& v: vs) {
        v /= 2.0;
        v += Vector3d{0.1, 0.24, 0.0};
    }

    Polygon quad(vs);

    TEST_ASSERT_CLOSE(quad.area(), 1.0);

    Vector3d cell_c = quad.center();

    // Тестовый набор объемных долей
    std::array test_alpha = {1.3e-9, 0.123, 0.253, 0.5, 0.78, 0.934, 1.0 - 2.3e-9};

    // Локальное число Куранта (тестовый набор)
    std::array test_CFL = {3.8e-3, 0.143, 0.254, 0.3634, 0.5, 0.6423};

    // Тестовый набор нормалей
    std::array test_n = {
        Vector3d{1.0, 0.0, 0.0},
        Vector3d{0.0, 1.0, 0.0},
        Vector3d{0.2, 0.4, 0.0},
        Vector3d{-0.2, 0.4, 0.0}
    };
    for (auto& n: test_n) {
        n.normalize();
    }

    for (double alpha: test_alpha) {
        for (double CFL: test_CFL) {
            for (const Vector3d& ni: test_n) {
                Vector3d p = cell_c + quad_find_section(alpha, ni) * ni;
                for (int i = 0; i < 4; ++i) {
                    Vector3d fc = 0.5 * (vs[i] + vs[(i + 1) % 4]);
                    Vector3d face_n = (fc - quad.center()).normalized();

                    Polygon clip = quad.clip(fc + 2 * CFL * (cell_c - fc), -face_n);

                    TEST_ASSERT_CLOSE(clip.area(), CFL);

                    double a_sig1 = clip.clip(p, ni).area() / clip.area();
                    double a_sig2 = average_flux(alpha, ni.dot(face_n), CFL);
                    double a_sig3 = average_flux(alpha, ni, face_n, CFL);

                    TEST_ASSERT_CLOSE(a_sig1, a_sig2);
                    TEST_ASSERT_EQ(a_sig2, a_sig3);
                }
            }
        }
    }
}

TEST_CASE("average_flux 3D") {
    Polyhedron cube = Polyhedron::Cube();
    cube.move({0.1, 0.23, 0.77});

    TEST_ASSERT_CLOSE(cube.volume(), 1.0);

    Vector3d cell_c = cube.center();

    // Тестовый набор объемных долей
    std::array test_alpha = {1.3e-9, 0.123, 0.253, 0.5, 0.78, 0.934, 1.0 - 2.3e-9};

    // Локальное число Куранта (тестовый набор)
    std::array test_CFL = {3.8e-3, 0.143, 0.254, 0.3634, 0.5, 0.6423};

    // Тестовый набор нормалей
    std::array test_n = {
        Vector3d{1.0, 0.0, 0.0},
        Vector3d{0.0, 1.0, 0.0},
        Vector3d{0.0, 0.0, 1.0},
        Vector3d{0.2, 0.4, 0.3},
        Vector3d{-0.2, 0.4, 0.9}
    };
    for (auto& n: test_n) {
        n.normalize();
    }

    for (double alpha: test_alpha) {
        for (double CFL: test_CFL) {
            for (const Vector3d& ni: test_n) {
                Vector3d p = cell_c + cube_find_section(alpha, ni) * ni;
                for (int i = 0; i < 6; ++i) {
                    Vector3d fc = cube.face_center(i);
                    Vector3d face_n = cube.face_normal(i);

                    Polyhedron clip = cube.clip(fc + 2 * CFL * (cell_c - fc), -face_n);

                    TEST_ASSERT_CLOSE(clip.volume(), CFL);

                    double a_sig1 = clip.clip(p, ni).volume() / clip.volume();
                    double a_sig2 = average_flux(alpha, ni, face_n, CFL);

                    TEST_ASSERT_CLOSE(a_sig1, a_sig2);
                }
            }
        }
    }
}

int main() {
    return TestRegistry::instance().run();
}