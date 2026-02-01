/// @file mmf_bubble_collapse.cpp
/// @brief Bubble collapse under the water
///
/// Nourgaliev R., Dinh T., Theofanous T. Adaptive characteristics-based matching
/// for compressible multifluid dynamics // Journal of Computational Physics.
/// –– 2006. –– Vol. 213, no. 2. –– P. 500–529.

#include <iomanip>

#include <zephyr/geom/primitives/polygon.h>
#include <zephyr/geom/generator/rectangle.h>
#include <zephyr/mesh/euler/eu_mesh.h>

#include <zephyr/phys/literals.h>
#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/matter/eos/stiffened_gas.h>
#include <zephyr/math/solver/mm_fluid.h>

#include <zephyr/io/pvd_file.h>
#include <zephyr/utils/mpi.h>
#include <zephyr/utils/threads.h>
#include <zephyr/utils/stopwatch.h>

using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::math::mmf;

using zephyr::io::PvdFile;
using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;

void init_cells(EuMesh& mesh, const MixturePT& mixture, Storable<PState> z) {
    PState z_air(
            1.0_kg_m3,          // density
            Vector3d::Zero(),   // velocity
            1.0_bar,            // pressure
            Fractions::Pure(0), // mass fractions
            mixture);

    PState z_water(
            1000.0_kg_m3,       // density
            Vector3d::Zero(),   // velocity
            1.0_bar,            // pressure
            Fractions::Pure(1), // mass fractions
            mixture);

    PState z_shock(
            1323.65_kg_m3,          // density
            {681.58_m_s, 0.0, 0.0}, // velocity
            1.9e4_bar,              // pressure
            Fractions::Pure(1),     // mass fractions
            mixture);

    double r = 3.0_mm;
    double x_bubble = 6.0_mm;
    double y_bubble = 6.0_mm;
    double x_shock  = 0.6_mm;

    Vector3d bubble_center = {x_bubble, y_bubble, 0};
    auto in_water = [bubble_center, r](const Vector3d &v) -> bool {
        return (v - bubble_center).norm() > r;
    };

    mesh.for_each([&](EuCell &cell) {
        if (cell.center().x() < x_shock) {
            cell[z] = z_shock;
            return;
        }

        double vol_frac1 = cell.approx_vol_fraction(in_water);
        if (vol_frac1 == 0.0 || vol_frac1 == 1.0) {
            // Pure cell (air or water)
            cell[z] = vol_frac1 > 0.5 ? z_water : z_air;
        }
        else {
            // Mixed cell (air and water)
            double vol_frac0 = cell.polygon().disk_clip_area(bubble_center, r) / cell.volume();
            vol_frac1 = 1.0 - vol_frac0;

            mmf::PState &z0 = z_air;
            mmf::PState &z1 = z_water;

            // rho = sum a_i rho_i
            double    density   = vol_frac0 * z0.density + vol_frac1 * z1.density;
            Fractions mass_frac = {vol_frac0 * z0.density / density, vol_frac1 * z1.density / density};
            mass_frac.normalize();
            Vector3d  velocity  = mass_frac[0] * z0.velocity + mass_frac[1] * z1.velocity;
            double    pressure  = vol_frac0 * z0.pressure + vol_frac1 * z1.pressure;

            PState mix(density, velocity, pressure, mass_frac, mixture);

            cell[z] = mix;
        }
    });
}

int main(int argc, char** argv) {
    mpi::handler handler(argc, argv);
    threads::init(argc, argv);
    threads::info();
    threads::off();

    // Generator of a Cartesian grid
    generator::Rectangle gen(0.0_cm, 1.2_cm, 0.0_cm, 1.2_cm);
    gen.set_nx(120);
    gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::WALL,
                        .bottom=Boundary::WALL, .top=Boundary::WALL});

    // Create mesh
    EuMesh mesh(gen);

    // Create EoS of materials and mixture
    Eos::Ptr air = IdealGas::create("Air");
    Eos::Ptr water = StiffenedGas::create("Water");
    MixturePT mixture = {air, water};

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
    mesh.set_decomposition("XY");
    mesh.set_max_level(0);
    mesh.set_distributor(solver.distributor());

    // Files for output
    PvdFile pvd("BC", "output");
    PvdFile pvd_bubble("bubble", "output");

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
        init_cells(mesh, mixture, data.init);
        solver.set_flags(mesh);
        mesh.refine();
    }
    init_cells(mesh, mixture, data.init);

    size_t n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;
    double max_time = 2.0 * 4.5_us;

    Stopwatch elapsed(true);
    while (curr_time < max_time) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << 1.0e6 * curr_time << " us\n";
            pvd.save(mesh, curr_time);

            solver.interface_recovery(mesh);
            auto bubble = solver.domain(mesh, 0);
            pvd_bubble.save(bubble, curr_time);

            next_write += 0.0; //max_time / 50;
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
    auto bubble = solver.domain(mesh, 0);
    pvd_bubble.save(bubble, max_time);

    elapsed.stop();

    std::cout << "\nElapsed:      " << elapsed.extended_time()
              << " ( " << elapsed.milliseconds() << " ms)\n";

    return 0;
}
