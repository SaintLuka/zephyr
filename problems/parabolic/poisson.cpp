#include <cmath>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/pvd_file.h>

using namespace zephyr::geom;
using namespace zephyr::mesh;
using namespace zephyr::io;
using generator::Rectangle;

double left_x(double y) {
    return 0.0;
}

double right_x(double y) {
    return std::sin(5) * std::cos(M_PI * std::pow(y, 2) / 2);
}

double left_y(double x) {
    return std::sin(5 * x);
}

double right_y(double x) {
    return - std::sin(5 * x);
}

double func(double x, double y) {
    return std::sin(5 * x) * (
        (25 + std::pow(M_PI * y, 2)) * std::cos(M_PI * std::pow(y, 2) / 2) +
        M_PI * std::sin(M_PI * std::pow(y, 2) / 2)
    );
}

void OpL(EuMesh &mesh, Storable<double> &in_var, Storable<double> &out_var) {
    for (auto cell: mesh) {
        cell[out_var] = 0.0;
    }
    for (auto cell: mesh) {
        if (cell.x() - cell.hx() / 2 < 1e-11) {
            cell[out_var] += 2 * cell[in_var] / std::pow(cell.hx(), 2);
        } else {
            cell[out_var] += (cell[in_var] - cell.face(Side2D::LEFT).neib(in_var)) / std::pow(cell.hx(), 2);
        }
        if (cell.x() + cell.hx() / 2 > 1 - 1e-11) {
            cell[out_var] += 2 * cell[in_var] / std::pow(cell.hx(), 2);
        } else {
            cell[out_var] += (cell[in_var] - cell.face(Side2D::RIGHT).neib(in_var)) / std::pow(cell.hx(), 2);
        }
        if (cell.y() - cell.hy() / 2 < 1e-11) {
            cell[out_var] += 2 * cell[in_var] / std::pow(cell.hy(), 2);
        } else {
            cell[out_var] += (cell[in_var] - cell.face(Side2D::BOTTOM).neib(in_var)) / std::pow(cell.hy(), 2);
        }
        if (cell.y() + cell.hy() / 2 > sqrt(2) - 1e-11) {
            cell[out_var] += 2 * cell[in_var] / std::pow(cell.hy(), 2);
        } else {
            cell[out_var] += (cell[in_var] - cell.face(Side2D::TOP).neib(in_var)) / std::pow(cell.hy(), 2);
        }
    }
}

void RHS(EuMesh &mesh, Storable<double> &rhs) {
    for (auto cell: mesh) {
        cell[rhs] = func(cell.x(), cell.y());
        if (cell.x() - cell.hx() / 2 < 1e-11) {
            cell[rhs] += 2 * left_x(cell.y()) / std::pow(cell.hx(), 2);
        }
        if (cell.x() + cell.hx() / 2 > 1 - 1e-11) {
            cell[rhs] += 2 * right_x(cell.y()) / std::pow(cell.hx(), 2);
        }
        if (cell.y() - cell.hy() / 2 < 1e-11) {
            cell[rhs] += 2 * left_y(cell.x()) / std::pow(cell.hy(), 2);
        }
        if (cell.y() + cell.hy() / 2 > sqrt(2) - 1e-11) {
            cell[rhs] += 2 * right_y(cell.x()) / std::pow(cell.hy(), 2);
        }
    }
}

double get_err(EuMesh &mesh, Storable<double> &lu_cur, Storable<double> &rhs) {
    double err = std::abs(mesh[0][lu_cur] - mesh[0][rhs]);
    for (auto cell: mesh) {
        double cur_err = std::abs(cell[lu_cur] - cell[rhs]);
        if (cur_err > err) {
            err = cur_err;
        }
    }
    return err;
}

double aim_func(double x, double y) {
    return std::sin(5 * x) * std::cos(M_PI * std::pow(y, 2) / 2);
}

void gd(EuMesh &mesh, Storable<double> &u_cur, double rtol, int max_iter, PvdFile &pvd) {
    for (auto cell: mesh) {
        cell[u_cur] = 0.0;
    }

    Storable<double> lu_cur = mesh.add<double>("lu_cur");
    Storable<double> rhs = mesh.add<double>("rhs");

    RHS(mesh, rhs);

    int counter = 0;
    double tau = std::pow(mesh[0].hx(), 2) / 16;

    do {
        OpL(mesh, u_cur, lu_cur);
        for (auto cell: mesh) {
            cell[u_cur] -= tau * (cell[lu_cur] - cell[rhs]);
        }
        counter++;
        if (counter % 1000 == 0) {
            std::cout << "Iteration: " << counter << std::endl;
        }
    } while ((get_err(mesh, lu_cur, rhs) > rtol) && (counter < max_iter));
    pvd.save(mesh, counter);
}

int main() {
    Rectangle gen(0.0, 1.0, 0.0, sqrt(2));  // [x_min, x_max, y_min, y_max]
    gen.set_sizes(200, 200);
    gen.set_boundaries({.left=Boundary::UNDEFINED, .right=Boundary::UNDEFINED,
                        .bottom=Boundary::UNDEFINED, .top=Boundary::UNDEFINED});
    EuMesh mesh(gen);

    Storable<double> u_cur = mesh.add<double>("u_cur");
    Storable<double> u_orig = mesh.add<double>("u_orig");

    PvdFile pvd("poisson", "output");
    pvd.variables.append("u_approx", u_cur);
    pvd.variables.append("u_orig", u_orig);

    for (auto cell: mesh) {
        cell[u_orig] = aim_func(cell.x(), cell.y());
    }

    gd(mesh, u_cur, 1e-10, 10000, pvd);
    return 0;
}