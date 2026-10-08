/// @file mmf_riemann_2D.cpp
/// @brief Two-dimensional Riemann problem (four materials).
/// Kurganov, Alexander, Tadmor, Eitan. Solution of two-dimensional Riemann
/// problems for gas dynamics without Riemann problem solvers.
/// Numer. Methods Partial Differ. Equ. 18(5), 2002. — c. 584-608.

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

#include <zephyr/phys/tests/test_2D.h>

using namespace zephyr::io;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::mmf;

using zephyr::mesh::EuMesh;

using zephyr::utils::Stopwatch;
using zephyr::utils::threads;
using zephyr::utils::mpi;
using zephyr::phys::Riemann2D;

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

    // Test problem
    Riemann2D test(5);

    // Generator of a Cartesian grid
    Rectangle gen(test.xmin(), test.xmax(), test.ymin(), test.ymax());
    gen.set_boundaries(test.boundaries());
    gen.set_nx(200);

    // Create mesh
    EuMesh mesh(gen);

    // Create EoS of materials and mixture
    MixturePT mixture;
    mixture += test.get_eos(0);
    mixture += test.get_eos(1);
    mixture += test.get_eos(2);
    mixture += test.get_eos(3);

    // Create and configure solver
    MmFluid solver(mixture);
    solver.set_CFL(0.5);
    solver.set_accuracy(2);
    solver.set_limiter("van Albada");
    solver.set_method(Fluxes::CRP);
    solver.set_crp_mode(CrpMode::PLIC);
    solver.set_splitting(DirSplit::SIMPLE);

    // Add data fields, choose main data layer
    auto data = solver.add_types(mesh);
    auto z = data.init;

    // Configure mesh
    mesh.set_max_level(3);
    mesh.set_decomposition("XY");
    mesh.set_distributor(solver.distributor());

    // Задание начальных данных
    auto init_cells = [&test, z, mixture](EuMesh& mesh) {
        mesh.for_each([&](EuCell& cell) {
            Vector3d r = cell.center();

            double    rho  = test.density(r);
            Vector3d  v    = test.velocity(r);
            double    P    = test.pressure(r);
            Fractions beta = test.fractions(r);

            cell[z] = PState(rho, v, P, beta, mixture);
        });
    };

    // Files for output
    PvdFile pvd("mesh", "output");
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

    // Initial conditions (adaptive to initial data)
    for (int k = 0; mesh.adaptive() && (k < mesh.max_level() + 3); ++k) {
        init_cells(mesh);
        solver.set_flags(mesh);
        mesh.refine();
    }
    init_cells(mesh);

    size_t n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;

    Stopwatch elapsed(true);
    while (curr_time < test.max_time()) {
        if (curr_time >= next_write) {
            mpi::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";
            pvd.save(mesh, curr_time);

            solver.interface_recovery(mesh);
            for (int i = 0; i < mixture.size(); ++i) {
                auto domain = solver.domain(mesh, i);
                pvd_domains[i].save(domain, curr_time);
            }

            next_write += test.max_time() / 50;
        }

        // Finish exactly at max_time
        solver.set_max_dt(test.max_time() - curr_time);

        // Integration step
        solver.update(mesh);
        solver.set_flags(mesh);
        mesh.refine();

        curr_time += solver.dt();
        n_step += 1;
    }
    pvd.save(mesh, test.max_time());

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
