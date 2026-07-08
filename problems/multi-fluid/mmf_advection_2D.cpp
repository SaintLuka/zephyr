/// @file mmf_advection_2D.cpp
/// @brief Two-dimensional advection problem for two materials.
/// Problem is simulated with multimaterial hydrodynamic solver.

#include <iomanip>

#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/matter/eos/stiffened_gas.h>
#include <zephyr/phys/matter/eos/mie_gruneisen.h>
#include <zephyr/phys/matter/mixture_pt.h>
#include <zephyr/math/solver/mm_fluid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/utils/threads.h>

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::mmf;

using zephyr::io::PvdFile;
using zephyr::utils::threads;

// Indicator function of the second material
bool inside_circle(const Vector3d &v) {
    Vector3d c = {0.15, 0.5, 0};
    double r = 0.1;
    return (v - c).norm() < r;
}

// Indicator function of the second material
bool inside_square(const Vector3d &v) {
    Vector3d c = {0.15, 0.5, 0};
    double l = 0.1;
    return (v - c).cwiseAbs().maxCoeff() < l;
}

// Choose function
auto inside = inside_circle;

int main() {
    threads::on();

    // Two identical/different materials
    Eos::Ptr eos1 = IdealGas::create("Air");
    //Eos::Ptr eos2 = IdealGas::create("Air");
    Eos::Ptr eos2 = StiffenedGas::create("Water");
    //Eos::Ptr eos2 = MieGruneisen::create("Fe");

    // Formal mixture
    MixturePT mixture = {eos1, eos2};

    // Create and configure solver
    MmFluid solver(mixture);
    solver.set_CFL(0.5);
    solver.set_accuracy(1);
    solver.set_method(Fluxes::CRP);
    solver.set_crp_mode(CrpMode::PLIC);
    solver.set_splitting(DirSplit::SIMPLE);

    // Generator of a Cartesian grid
    generator::Rectangle gen(0.0, 1.0, 0.0, 0.7);
    gen.set_nx(200);
    gen.set_boundaries({Boundary::ZOE, Boundary::ZOE, Boundary::ZOE, Boundary::ZOE});

    // Create mesh
    EuMesh mesh(gen);

    // Add data fields, choose main data layer
    auto data = solver.add_types(mesh);
    auto z = data.init;

    // Files for output
    PvdFile pvd("mesh", "output");
    PvdFile pvd_body("body", "output");

    // Variables to save
    pvd.variables = {"level"};
    pvd.variables += {"cln", [z](EuCell cell) -> double { return cell[z].mass_frac.index(); }};
    pvd.variables += {"rho", [z](EuCell cell) -> double { return cell[z].density; }};
    pvd.variables += {"vx",  [z](EuCell cell) -> double { return cell[z].velocity.x(); }};
    pvd.variables += {"vy",  [z](EuCell cell) -> double { return cell[z].velocity.y(); }};
    pvd.variables += {"e",   [z](EuCell cell) -> double { return cell[z].energy; }};
    pvd.variables += {"P",   [z](EuCell cell) -> double { return cell[z].pressure; }};
    pvd.variables += {"T",   [z](EuCell cell) -> double { return cell[z].temperature; }};
    pvd.variables += {"b0",  [z](EuCell cell) -> double { return cell[z].mass_frac[0]; }};
    pvd.variables += {"b1",  [z](EuCell cell) -> double { return cell[z].mass_frac[1]; }};
    pvd.variables += {"a0",  [z](EuCell cell) -> double { return cell[z].alpha(0); }};
    pvd.variables += {"a1",  [z](EuCell cell) -> double { return cell[z].alpha(1); }};
    pvd.variables += {"rho0",[z](EuCell cell) -> double { return cell[z].densities[0]; }};
    pvd.variables += {"rho1",[z](EuCell cell) -> double { return cell[z].densities[1]; }};
    pvd.variables += {"n.x", [n=data.n](EuCell cell) -> double { return cell[n][0].x(); }};
    pvd.variables += {"n.y", [n=data.n](EuCell cell) -> double { return cell[n][0].y(); }};

    mesh.for_each([&](EuCell &cell) {
        const Vector3d V0 = {70.0, -35.0, 0.0};

        const PState z1(
                eos1->density(),    // density
                V0,                 // velocity
                1.0e5,              // pressure
                Fractions::Pure(0), // mass fractions
                mixture);

        const PState z2(
                eos2->density(),    // density
                V0,                 // velocity
                1.0e5,              // pressure
                Fractions::Pure(1), // mass fractions
                mixture);

        double vol_frac2 = cell.approx_vol_fraction(inside);
        if (vol_frac2 == 0.0 || vol_frac2 == 1.0) {
            // Pure cell
            cell[z] = vol_frac2 < 0.5 ? z1 : z2;
        }
        else {
            // Mixed cell
            vol_frac2 = cell.volume_fraction(inside, 10'000);
            double vol_frac1 = 1.0 - vol_frac2;

            cell[z] = PState::Mix1(mixture, {vol_frac1, vol_frac2}, {z1, z2});
        }
    });

    size_t n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;
    double max_time = 0.01;

    while (curr_time < max_time) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);

            solver.interface_recovery(mesh);
            auto body = solver.domain(mesh, 1);
            pvd_body.save(body, curr_time);

            next_write += max_time / 100;
        }

        // Finish exactly at max_time
        solver.set_max_dt(max_time - curr_time);

        // Integration step
        solver.update(mesh);

        curr_time += solver.dt();
        n_step += 1;
    }
    pvd.save(mesh, max_time);

    solver.interface_recovery(mesh);
    auto body = solver.domain(mesh, 1);
    pvd_body.save(body, curr_time);

    return 0;
}
