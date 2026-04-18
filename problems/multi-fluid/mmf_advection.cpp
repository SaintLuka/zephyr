/// @file mmf_advection.cpp
/// @brief Two-dimensional advection problem for two materials.
/// Problem is simulated with multimaterial hydrodynamic solver.

#include <iomanip>

#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/matter/mixture_pt.h>
#include <zephyr/math/solver/mm_fluid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/utils/threads.h>

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::mmf;

using zephyr::io::PvdFile;
using zephyr::utils::threads;

// Indicator function of the domain
bool inside(const Vector3d& r) {
    return std::abs(r.x() - 0.15) < 0.1 && std::abs(r.y() - 0.50) < 0.1;
}

int main() {
    threads::on();

    // Two identical materials
    double gamma = 1.4;
    Eos::Ptr sg1 = IdealGas::create(gamma, 1.0);
    Eos::Ptr sg2 = IdealGas::create(gamma, 1.0);

    // Formal mixture
    MixturePT mixture = {sg1, sg2};

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
    PvdFile pvd("Advection2D", "output");
    PvdFile pvd_domain("domain", "output");

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

    // Initial conditions
    for (auto cell: mesh) {
        cell[z].velocity    = {0.7, -0.35, 0.0};
        cell[z].density     = 1.0;
        cell[z].pressure    = 1.0 / gamma;
        //cell[z].energy      = 1.0 / (gamma * (gamma - 1.0));
        //cell[z].temperature = 1.0;

        //if (r.x() < 0.1) cell[z].mass_frac[0]  = 0.0;
        //if (std::abs(r.x() - 0.2) < 0.1) cell[z].mass_frac[0]  = (r.x() - 0.1) / 0.2;
        //if (r.x() > 0.3) cell[z].mass_frac[0]  = 1.0;

        cell[z].mass_frac[0] = inside(cell.center()) ? 1.0 : 0.0;
        cell[z].mass_frac[1] = 1.0 - cell[z].mass_frac[0];

        cell[z].densities[0] = cell[z].mass_frac[0] > 0.0 ? cell[z].density : NAN;
        cell[z].densities[1] = cell[z].mass_frac[1] > 0.0 ? cell[z].density : NAN;

        //cell[z].density = 1.0 / mixture.volume_PT(cell[z].pressure, cell[z].temperature, cell[z].mass_frac);
        //cell[z].energy = mixture.energy_PT(cell[z].pressure, cell[z].temperature, cell[z].mass_frac);

        cell[z].energy      = mixture.energy_rP     (cell[z].density, cell[z].pressure, cell[z].mass_frac);
        cell[z].temperature = mixture.temperature_rP(cell[z].density, cell[z].pressure, cell[z].mass_frac);
    }

    size_t n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;
    double max_time = 1.0;

    while (curr_time < max_time) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);

            solver.interface_recovery(mesh);
            auto domain = solver.domain(mesh, 0);
            pvd_domain.save(domain, curr_time);

            next_write += max_time / 200;
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
    auto domain = solver.domain(mesh, 0);
    pvd_domain.save(domain, curr_time);

    return 0;
}
