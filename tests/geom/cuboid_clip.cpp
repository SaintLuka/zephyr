#include <iostream>
#include <zephyr/geom/geom.h>
#include <zephyr/geom/sections.h>
#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/math/calc/roots.h>
#include <zephyr/utils/numpy.h>
#include <zephyr/utils/pyplot.h>

using namespace zephyr;
using namespace zephyr::geom;
using namespace zephyr::mesh;

// Размеры параллелепипеда
const double A = 1.2;
const double B = 1.4;
const double C = 1.0;

void test1(Polyhedron poly, Vector3d normal) {
    double p_min1 = +1.0e100;
    double p_max1 = -1.0e100;
    for (int i = 0; i < poly.n_verts(); ++i) {
        p_min1 = std::min(p_min1, poly.vertex(i).dot(normal));
        p_max1 = std::max(p_max1, poly.vertex(i).dot(normal));
    }
    double p_min = p_min1 - 0.05 * (p_max1 - p_min1);
    double p_max = p_max1 + 0.05 * (p_max1 - p_min1);

    int N = 1001;

    std::vector<double> ps(N);
    std::vector<double> Vs1(N);
    std::vector<double> Vs2(N);
    std::vector<double> err(N);
    for (int i = 0; i < N; ++i) {
        Vector3d p = (p_min + (p_max - p_min) * i / N) * normal;

        ps[i] = p.dot(normal);

        Polyhedron clip = poly.clip(p, normal);
        Vs1[i] = clip.volume() / poly.volume();
        Vs2[i] = cube_volume_fraction(normal, p.dot(normal), A, B, C);

        err[i] = std::abs(Vs1[i] - Vs2[i]);
    }

    utils::pyplot plt;
    plt.figure({.figsize={9.0, 4.0}, .dpi=170});

    plt.subplot(1, 2, 1);
    plt.title("$\\alpha(p)$");
    plt.grid(true);
    plt.plot(ps, Vs1, {.linestyle="solid", .color="blue", .label="clip_volume"});
    plt.plot(ps, Vs2, {.linestyle="dotted", .color="orange", .label="formula"});

    plt.subplot(1, 2, 2);
    plt.title("$Error$");
    plt.grid(true);
    plt.plot(ps, err, "k");

    plt.tight_layout();
    plt.show();
}

void test2(Polyhedron poly, Vector3d normal) {
    int N = 501;
    auto alpha = np::linspace(0.0, 1.0, N);

    auto ps  = np::zeros(N);
    auto err = np::zeros(N);
    for (int i = 0; i < N; ++i) {
        ps[i]  = cube_find_section(normal, alpha[i], A, B, C);
        err[i] = std::abs(alpha[i] - cube_volume_fraction(normal, ps[i], A, B, C));
    }

    utils::pyplot plt;
    plt.figure({.figsize={9.0, 4.0}, .dpi=170});

    plt.subplot(1, 2, 1);
    plt.title("$p(\\alpha)$");
    plt.grid(true);
    plt.plot(alpha, ps, {.color="blue", .label="clip_volume"});

    plt.subplot(1, 2, 2);
    plt.title("$Error$");
    plt.grid(true);
    plt.plot(alpha, err, "k");

    plt.tight_layout();
    plt.show();
}

int main() {
    auto poly = Polyhedron::Cuboid(A, B, C);

    Vector3d n = {-0.12, 0.16, -0.51};
    //Vector3d n = {0.51, 0.53, 0.52};
    n.normalize();

    test1(poly, n);
    test2(poly, n);

    return 0;
}