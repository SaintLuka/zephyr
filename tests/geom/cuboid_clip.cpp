#include <iostream>
#include <zephyr/geom/geom.h>
#include <zephyr/geom/sections.h>

using namespace zephyr::geom;

int main() {
    Vector3d n = {0.1, 0.2, 0.3};
    n.normalize();

    double p = 0.1;
    Vector3d P = p * n;

    auto poly = Polyhedron::Cuboid(1.0, 1.0, 1.0);

    double V1 = poly.clip(P, n).volume();
    double V2 = cube_volume_fraction(n, p, 1.0, 1.0, 1.0);

    std::cout << V1 << " " << V2 << "\n";

    return 0;
}