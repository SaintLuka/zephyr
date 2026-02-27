/// @file advection.cpp
/// @brief Numerical solution of a scalar advection problem.
/// A simple upwind scheme of the first order accuracy is used.

#include <iomanip>

#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/pvd_file.h>

using zephyr::geom::Vector3d;
using zephyr::geom::Box;
using zephyr::geom::Boundary;
using zephyr::geom::generator::Rectangle;
using zephyr::mesh::Storable;
using zephyr::mesh::EuMesh;
using zephyr::mesh::EuCell;
using zephyr::io::PvdFile;

// Velocity vector field
Vector3d velocity(const Vector3d& c) {
    return { 1.0, 0.3 + 0.3*std::sin(4 * M_PI * c.x()), 0.0 };
}

int main() {
    // Generator of a Cartesian grid
    Rectangle rect(0.0, 1.0, 0.0, 0.6, true);
    rect.set_nx(200);
    rect.set_boundaries({
        .left   = Boundary::PERIODIC, .right = Boundary::PERIODIC,
        .bottom = Boundary::PERIODIC, .top   = Boundary::PERIODIC});

    // Create mesh
    EuMesh mesh(rect);

    // Add data fields
    auto u1 = mesh.add<double>("u1");
    auto u2 = mesh.add<double>("u2");

    // Files for output
    PvdFile pvd("mesh", "output");

    // Variables to save
    pvd.variables.append("u", u1);
    pvd.variables += {"vx", [](EuCell cell) -> double { return velocity(cell.center()).x(); } };
    pvd.variables += {"vy", [](EuCell cell) -> double { return velocity(cell.center()).y(); } };

    // Initial conditions
    Box box = mesh.bbox();
    Vector3d vc = box.center();
    double D = 0.1 * box.diameter();
    for (auto cell: mesh) {
        cell[u1] = (cell.center() - vc).norm() < D ? 1.0 : 0.0;
        cell[u2] = 0.0;
    }

    // Courant number (< 1.0)
    double CFL = 0.5;

    int n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    while(curr_time <= 1.0) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);
            next_write += 0.02;
        }

        // Estimate timestep
        double dt = std::numeric_limits<double>::max();
        for (auto cell: mesh) {
            double max_area = 0.0;
            for (auto &face: cell.faces()) {
                max_area = std::max(max_area, face.area());
            }
            double dx = cell.volume() / max_area;
            dt = std::min(dt, dx / velocity(cell.center()).norm());
        }
        dt *= CFL;

        // Integrate using upwind scheme
        for (auto cell: mesh) {
            double zc = cell[u1];

            double fluxes = 0.0;
            for (auto& face: cell.faces()) {
                double zn = face.neib(u1);

                double af = velocity(face.center()).dot(face.normal());
                double a_p = std::max(af, 0.0);
                double a_m = std::min(af, 0.0);

                fluxes += (a_p * zc + a_m * zn) * face.area();
            }

            cell[u2] = zc - dt * fluxes / cell.volume();
        }

        // Update time layers
        for (auto cell: mesh) {
            cell[u1] = cell[u2];
            cell[u2] = 0.0;
        }

        n_step += 1;
        curr_time += dt;
    }

    return 0;
}