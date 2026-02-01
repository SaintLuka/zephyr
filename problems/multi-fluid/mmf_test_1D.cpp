/// @file mmf_test_1D.cpp
/// @brief One-dimensional Riemann problems for multimaterial hydrodynamic.

#include <iostream>
#include <iomanip>

#include <zephyr/geom/generator/strip.h>
#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/phys/tests/test_1D.h>
#include <zephyr/math/solver/mm_fluid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/io/csv_file.h>

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::mmf;

using zephyr::io::PvdFile;
using zephyr::io::CsvFile;
using zephyr::utils::threads;

int main() {
    // Test problems
    //RarefiedWater test;
    //Multimat1D test(1, 1, 0);
    ToroTest test(1, true);
    //test.adjust_cv();
    
    // Generator of quasi one-dimensional mesh
    generator::Strip gen(test.xmin(), test.xmax());
    gen.set_size(1000);

    // Create mesh
    EuMesh mesh(gen);

    // Test class provides materials
    MixturePT mixture = test.mixture_PT();

    // Create and configure solver
    MmFluid solver(mixture);
    solver.set_CFL(0.5);
    solver.set_accuracy(2);
    solver.set_method(Fluxes::HLLC);
    
    // Add data fields, choose main data layer
    auto data = solver.add_types(mesh);
    auto z = data.init;
    
    double curr_time = 0.0;

    // Files for output
    PvdFile pvd("test1D", "output");

    // Variables to save
    pvd.variables += {"cln", [z](EuCell cell) -> double { return cell[z].mass_frac.index(); }};
    pvd.variables += {"rho", [z](EuCell cell) -> double { return cell[z].density; }};
    pvd.variables += {"u",   [z](EuCell cell) -> double { return cell[z].velocity.x(); }};
    pvd.variables += {"e",   [z](EuCell cell) -> double { return cell[z].energy; }};
    pvd.variables += {"P",   [z](EuCell cell) -> double { return cell[z].pressure; }};
    pvd.variables += {"T",   [z](EuCell cell) -> double { return cell[z].temperature; }};
    pvd.variables += {"b0",  [z](EuCell cell) -> double { return cell[z].mass_frac[0]; }};
    pvd.variables += {"a0",  [z](EuCell cell) -> double { return cell[z].alpha(0); }};

    // Variables to save (exact solution)
    pvd.variables += {"rho.exact",
                      [&test, &curr_time](EuCell cell) -> double {
                          return test.density_t(cell.center(), curr_time);
                      }};
    pvd.variables += {"u.exact",
                      [&test, &curr_time](EuCell cell) -> double {
                          return test.velocity_t(cell.center(), curr_time).x();
                      }};
    pvd.variables += {"P.exact",
                      [&test, &curr_time](EuCell cell) -> double {
                          return test.pressure_t(cell.center(), curr_time);
                      }};
    pvd.variables += {"e.exact",
                      [&test, &curr_time](EuCell cell) -> double {
                          return test.energy_t(cell.center(), curr_time);
                      }};
    pvd.variables += {"T.exact",
                      [&test, &curr_time](EuCell cell) -> double {
                          return test.temperature_t(cell.center(), curr_time);
                      }};
    pvd.variables += {"b0.exact",
                      [&test, &curr_time](EuCell cell) -> double {
                          return test.fractions_t(cell.center(), curr_time)[0];
                      }};

    // Initial conditions
    for (auto cell: mesh) {
        Vector3d r = cell.center();
        cell[z].density  = test.density(r);
        cell[z].velocity = test.velocity(r);
        cell[z].pressure = test.pressure(r);
        cell[z].energy   = test.energy(r);

        cell[z].mass_frac = test.fractions(r);

        cell[z].densities[0] = test.fractions(r)[0] > 0.0 ? test.density(r) : NAN;
        cell[z].densities[1] = test.fractions(r)[1] > 0.0 ? test.density(r) : NAN;

        cell[z].temperature = mixture.temperature_rP(
                cell[z].density, cell[z].pressure, cell[z].mass_frac);
    }

    size_t n_step = 0;
    double next_write = 0.0;

    while (curr_time < test.max_time()) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);
            next_write += test.max_time() / 200;
        }

        // Finish exactly at max_time
        solver.set_max_dt(test.max_time() - curr_time);

        // Integration step
        solver.update(mesh);

        curr_time += solver.dt();
        n_step += 1;
    }

    // Save as csv datasets
    CsvFile csv("test1D.csv", 8, pvd.variables);
    csv.save(mesh);

    // Compute errors
    double r_err = 0.0, u_err = 0.0, p_err = 0.0, e_err = 0.0;
    double r_avg = 0.0, u_avg = 0.0, p_avg = 0.0, e_avg = 0.0;
    for (auto cell: mesh) {
        Vector3d r = cell.center();
        double V = cell.volume();

        r_err += V * std::abs(cell[z].density      - test.density_t(r, curr_time));
        u_err += V * std::abs(cell[z].velocity.x() - test.velocity_t(r, curr_time).x());
        p_err += V * std::abs(cell[z].pressure     - test.pressure_t(r, curr_time));
        e_err += V * std::abs(cell[z].energy       - test.energy_t(r, curr_time));

        r_avg += V * std::abs(cell[z].density     );
        u_avg += V * std::abs(cell[z].velocity.x());
        p_avg += V * std::abs(cell[z].pressure    );
        e_avg += V * std::abs(cell[z].energy      );
    }
    r_err /= r_avg;
    u_err /= u_avg;
    p_err /= p_avg;
    e_err /= e_avg;

    std::cout << "\nMean errors\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "    Density:  " << r_err << "\n";
    std::cout << "    Velocity: " << u_err << "\n";
    std::cout << "    Pressure: " << p_err << "\n";
    std::cout << "    Energy:   " << e_err << "\n";

    return 0;
}
