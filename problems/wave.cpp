/// @file wave.cpp
/// @brief Numerical solution of a wave equation.

#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>
#include <zephyr/io/pvd_file.h>

using namespace zephyr::geom;
using namespace zephyr::mesh;
using namespace zephyr::io;
using generator::Rectangle;

int main() {
    // Generator of a Cartesian grid
    Rectangle gen(-2.0, 1.0, -2.0, 1.0);
    gen.set_nx(400);
    gen.set_boundaries({.left=Boundary::WALL, .right=Boundary::WALL,
                        .bottom=Boundary::WALL, .top=Boundary::WALL});

    // Create mesh
    EuMesh mesh(gen);

    // Add data fields
    Storable<double> u_prev = mesh.add<double>("u_prev");
    Storable<double> u_curr = mesh.add<double>("u_curr");
    Storable<double> u_next = mesh.add<double>("u_next");

    // Files for output
    PvdFile pvd("wave", "output");
    pvd.variables.append("u", u_curr);

    // Initial conditions, assuming zero velocity
    for (auto cell: mesh) {
        double r = cell.center().norm();
        cell[u_prev] = 0.5 + 0.5 * std::tanh(50*(0.5 - r));
        cell[u_curr] = cell[u_prev];
    }

    // Estimate timestep
    double speed = 1.0;
    double dt = std::numeric_limits<double>::max();
    for (auto cell: mesh) {
        dt = std::min(dt, 0.5 * cell.incircle_diameter() / speed);
    }
    dt *= 0.95; // CFL condition

    // Main loop
    size_t n_step = 0;
    double curr_time = dt;
    while (curr_time < 5.0 && n_step < 5000) {
        // Save data to file
        if (n_step % 10 == 0) {
            std::cout << "Step: " << n_step << ";\tTime: " << curr_time << "\n";
            pvd.save(mesh, curr_time);
        }

        // Loop over the cells, compute new time layer
        for (auto cell: mesh) {
            double uc = cell[u_curr];

            double fluxes = 0.0;
            for (auto face: cell.faces()) {
                double un = face.neib(u_curr);
                if (face.is_boundary()) {
                    un = -uc;
                }

                double dr = (face.center() - cell.center()).norm();
                fluxes += (un - uc) * face.area() / (2.0 * dr);
            }

            cell[u_next] = 2.0 * cell[u_curr] - cell[u_prev] +
                std::pow(speed * dt, 2) * fluxes / cell.volume();
        }

        // Swap time layers
        for (auto cell: mesh) {
            cell[u_prev] = cell[u_curr];
            cell[u_curr] = cell[u_next];
        }

        // Update time and step
        curr_time += dt;
        ++n_step;
    }
    pvd.save(mesh, curr_time);
    return 0;
}