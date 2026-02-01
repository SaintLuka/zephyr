/// @file mmf_triple_point.cpp
/// @brief Triple point problem. Two or three material problem.
/// Shock wave in gases and Kelvin-Helmholtz instability.
///
/// [1] Pan S., Han L., Hu X., Adams N. A conservative interface-interaction
/// method for compressible multi-material flows // Journal of Computational
/// Physics. –– 2018. –– Vol. 371. –– P. 870–895.
/// [2] Dobrev V., Ellis T., Kolev T., Rieben R. High-order curvilinear finite
/// elements for axisymmetric lagrangian hydrodynamics // Computers and Fluids.
/// –– 2013. –– Vol. 83. –– P. 58–69.
/// [3] Galera S., Maire P.-H., Breil J. A two-dimensional unstructured cell-centered
/// multi-material ale scheme using vof interface reconstruction // Journal of
/// Computational Physics. –- 2010. –- Vol. 229, no. 16. –- P. 5755–5787.

#include <iomanip>

#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/matter/mixture_pt.h>
#include <zephyr/math/solver/mm_fluid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/utils/mpi.h>
#include <zephyr/utils/threads.h>
#include <zephyr/utils/stopwatch.h>

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::mmf;

using zephyr::mesh::EuMesh;

using zephyr::utils::Stopwatch;
using zephyr::utils::threads;
using zephyr::utils::mpi;

void init_cells(EuMesh& mesh, const MixturePT& mixture, Storable<PState> init) {
    double R_max = 1.0;
    double R_min = 0.125;
    double P_max = 1.0;
    double P_min = 0.1;

    double x_barrier = 1.0;
    double y_barrier = 1.5;

    int n = mixture.size() - 1;

    mesh.for_each([&](EuCell cell) {
        Vector3d v = cell.center();

        PState z = PState::Zero();
        if (v.x() > x_barrier && v.y() < y_barrier) {
            z.density   = R_max;
            z.pressure  = P_min;
            z.mass_frac = Fractions::Pure(0);
            z.densities = ScalarSet::PureNaN(0, R_max);
        }
        else if (v.x() < x_barrier) {
            z.density   = R_max;
            z.pressure  = P_max;
            z.mass_frac = Fractions::Pure(1);
            z.densities = ScalarSet::PureNaN(1, R_max);
        }
        else {
            // second or third material
            z.density   = R_min;
            z.pressure  = P_min;
            z.mass_frac = Fractions::Pure(n);
            z.densities = ScalarSet::PureNaN(n, R_min);
        }

        z.e() = mixture.energy_rP(z.rho(), z.P(), z.beta());
        z.T() = mixture.temperature_rP(z.rho(), z.P(), z.beta());

        cell[init] = z;
    });
}

int main(int argc, char** argv) {
    mpi::handler handler(argc, argv);
    threads::init(argc, argv);
    threads::info();
    //threads::off();

    // Generator of a Cartesian grid
    generator::Rectangle gen(0.0, 7.0, 0.0, 3.0);
    gen.set_nx(350);
    gen.set_boundaries({.left=Boundary::WALL, .right=Boundary::WALL,
                        .bottom=Boundary::WALL, .top=Boundary::WALL});

    // Create mesh
    EuMesh mesh(gen);

    // Create EoS of materials and mixture
    MixturePT mixture;
    mixture += IdealGas::create(1.4, 1.0);
    mixture += IdealGas::create(1.5, 1.0);
    // optional, same as the second material
    mixture += IdealGas::create(1.5, 1.0);

    // Create and configure solver
    MmFluid solver(mixture);
    solver.set_CFL(0.5);
    solver.set_accuracy(1);
    solver.set_method(Fluxes::CRP);
    solver.set_crp_mode(CrpMode::PLIC);
    solver.set_splitting(DirSplit::SIMPLE);

    // Add data fields, choose main data layer
    auto data = solver.add_types(mesh);
    auto z = data.init;

    // Configure mesh
    mesh.set_max_level(1);
    mesh.set_decomposition("XY");
    mesh.set_distributor(solver.distributor());

    // Files for output
    PvdFile pvd("TP", "output");
    std::vector<PvdFile> pvd_domains;
    for (int i = 0; i < mixture.size(); ++i) {
        pvd_domains.emplace_back("domain" + std::to_string(i), "output");
    }

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

    // Initial conditions (adaptive to initial data)
    for (int k = 0; mesh.adaptive() && k < mesh.max_level() + 3; ++k) {
        init_cells(mesh, mixture, z);
        solver.set_flags(mesh);
        mesh.refine();
    }
    init_cells(mesh, mixture, z);

    size_t n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;
    double max_time = 5.0;

    Stopwatch elapsed(true);
    while (curr_time < max_time) {
        if (curr_time >= next_write) {
            mpi::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);

            solver.interface_recovery(mesh);
            for (int i = 0; i < mixture.size(); ++i) {
                auto domain = solver.domain(mesh, i);
                pvd_domains[i].save(domain, curr_time);
            }

            next_write += max_time / 50;
        }

        // Finish exactly at max_time
        solver.set_max_dt(max_time - curr_time);

        // Integration step
        solver.update(mesh);
        solver.set_flags(mesh);
        mesh.refine();

        curr_time += solver.dt();
        n_step += 1;
    }
    pvd.save(mesh, max_time);

    solver.interface_recovery(mesh);
    for (int i = 0; i < mixture.size(); ++i) {
        auto domain = solver.domain(mesh, i);
        pvd_domains[i].save(domain, curr_time);
    }

    elapsed.stop();

    mpi::cout << "\nElapsed:      " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    return 0;
}
