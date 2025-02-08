/// @brief Тестирование многогранников, сечения, объемы и прочее

#include <zephyr/geom/polyhedron.h>
#include <zephyr/io/vtu_file.h>


using namespace zephyr::geom;
using namespace zephyr::io;


Polyhedron tetra() {
    std::vector<Vector3d> vs = {
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.3, 0.3, 0.5},
    };
    return Polyhedron(CellType::TETRA, vs);
}

Polyhedron pyramid() {
    std::vector<Vector3d> vs = {
            {0.0, 0.0, 0.0},
            {1.1, 0.0, 0.0},
            {0.9, 0.8, 0.0},
            {0.0, 1.0, 0.0},
            {0.6, 0.4, 0.8},
    };
    return Polyhedron(CellType::PYRAMID, vs);
}

Polyhedron wedge() {
    std::vector<Vector3d> vs = {
            {0.0, 0.0, 0.0},
            {0.9, 0.8, 0.0},
            {1.1, 0.0, 0.0},
            {0.1, 0.1, 0.5},
            {1.0, 0.9, 0.5},
            {1.2, 0.1, 0.5}
    };
    return Polyhedron(CellType::WEDGE, vs);
}

Polyhedron unit_cube_1() {
    std::vector<Vector3d> vs = {
            {-0.5, -0.5, -0.5},
            {+0.5, -0.5, -0.5},
            {+0.5, +0.5, -0.5},
            {-0.5, +0.5, -0.5},
            {-0.5, -0.5, +0.5},
            {+0.5, -0.5, +0.5},
            {+0.5, +0.5, +0.5},
            {-0.5, +0.5, +0.5},
    };

    Polyhedron poly(CellType::HEXAHEDRON, vs);

    return poly;
}

Polyhedron cuboid() {
    double a = 1.222;
    double b = 1.454;
    double c = 0.982;

    std::vector<Vector3d> vs = {
            {-0.5, -0.5, -0.5},
            {+0.5, -0.5, -0.5},
            {+0.5, +0.5, -0.5},
            {-0.5, +0.5, -0.5},
            {-0.5, -0.5, +0.5},
            {+0.5, -0.5, +0.5},
            {+0.5, +0.5, +0.5},
            {-0.5, +0.5, +0.5},
    };

    Polyhedron poly(CellType::HEXAHEDRON, vs);

    return poly;
}

double volfrac(double Pz, Vector3d n, double a, double b, double c) {
    std::vector<double> kek={std::abs(a * n.x()), std::abs(b * n.y()), std::abs(c * n.z())};
    std::sort(kek.begin(), kek.end());

    double xi = kek[0];
    double eta = kek[1];
    double chi = kek[2];

    std::cout << "  n: " << n << "\n";
    std::cout << "  xi : " << xi << "\n";
    std::cout << "  eta: " << eta << "\n";
    std::cout << "  chi: " << chi << "\n";

    double p = Pz + 0.5 * (xi + eta + chi);

    bool inv = p > 0.5 * (xi + eta + chi);

    if (inv) {
        std::cout << "  Inversed case\n";
        p = xi + eta + chi - p;
    }

    double vol = NAN;
    if (p < xi) {
        std::cout << "  case 1\n";
        vol = std::pow(p, 3) / (6.0 * xi * eta * chi);
    }
    else if (p < eta) {
        std::cout << "  case 2\n";
        vol = (3.0 * p * (p - xi) + xi * xi) / (6.0 * eta * chi);
    }
    else if (p < std::min(xi + eta, chi)) {
        std::cout << "  case 3\n";

        vol = (std::pow(p, 3) - std::pow(p - xi, 3) - std::pow(p - eta, 3)) / (6.0 * xi * eta * chi);
    }
    else {
        if (xi + eta < chi) {
            std::cout << "  case 4.1\n";
            vol = (2.0 * p - xi - eta) / (2.0 * chi);
        } else {
            std::cout << "  case 4.2\n";

            vol = (std::pow(p, 3) - std::pow(p - xi, 3) - std::pow(p - eta, 3) - std::pow(p - chi, 3)) / (6.0 * xi * eta * chi);
        }
    }

    if (inv) {
        vol = 1.0 - vol;
    }

    return vol;

}

int main() {

    AmrStorage cell(1);

    Polyhedron poly = unit_cube_1();
    //Polyhedron poly = unit_cube_1();

    double V =  poly.volume();

    std::cout << "Volume: " << V << "\n";

    Vector3d p = {0.2, 0.2, 0.2};
    Vector3d n = {0.1, 0.3, 0.5};
    n.normalize();

    cell[0] = AmrCell(poly);

    VtuFile::save("cell1.vtu", cell, Variables{});

    poly = poly.clip(p, n);

    cell[0] = AmrCell(poly);

    VtuFile::save("cell2.vtu", cell, Variables{});


    return 0;
}