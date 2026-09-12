/// @file mmf_advection_3D.cpp
/// @brief Three-dimensional advection problem for two materials.
/// Problem is simulated with multimaterial hydrodynamic solver.

#include <iomanip>

#include <zephyr/geom/generator/cuboid.h>
#include <zephyr/mesh/mesh.h>

#include <zephyr/phys/matter/eos/ideal_gas.h>
#include <zephyr/phys/matter/eos/stiffened_gas.h>
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
bool inside_ball(const Vector3d &v) {
    Vector3d c = Vector3d::Zero();
    double r = 0.5;
    return (v - c).norm() < r;
}

// Indicator function of the second material
bool inside_cube(const Vector3d &v) {
    Vector3d c = Vector3d::Zero();
    double l = 0.5;
    return (v - c).cwiseAbs().maxCoeff() < l;
}

// Choose function
auto inside = inside_ball;

// Setup initial conditions
void init_cells(Mesh &mesh, const MixturePT& mixture, Storable<PState> z) {
    auto eos1 = mixture[0];
    auto eos2 = mixture[1];

    mesh.for_each([&](Cell &cell) {
        const Vector3d V0 = {10.0, 5.0, 0.0};

        const PState z1(
                eos1.density(),     // density
                V0,                 // velocity
                1.0e5,              // pressure
                Fractions::Pure(0), // mass fractions
                mixture);

        const PState z2(
                eos2.density(),     // density
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
}

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
    generator::Cuboid gen(-2.0, 4.0, -2.0, 4.0, -2.0, 4.0);
    gen.set_nx(15);
    gen.set_boundaries({Boundary::ZOE, Boundary::ZOE, Boundary::ZOE,
                        Boundary::ZOE, Boundary::ZOE, Boundary::ZOE});

    // Create mesh
    Mesh mesh(gen);

    // Add data fields, choose main data layer
    auto data = solver.add_types(mesh);
    auto z = data.init;

    // Configure mesh
    mesh.set_decomposition("XY");
    mesh.set_max_level(2);
    mesh.set_distributor(solver.distributor());

    // Files for output
    PvdFile pvd("mesh", "output");
    PvdFile pvd_body("body", "output");
    pvd_body.options.polyhedral = true;

    // Variables to save
    pvd.variables = {"level"};
    pvd.variables += {"cln", [z](Cell cell) -> double { return cell[z].mass_frac.index(); }};
    pvd.variables += {"rho", [z](Cell cell) -> double { return cell[z].density; }};
    pvd.variables += {"vx",  [z](Cell cell) -> double { return cell[z].velocity.x(); }};
    pvd.variables += {"vy",  [z](Cell cell) -> double { return cell[z].velocity.y(); }};
    pvd.variables += {"e",   [z](Cell cell) -> double { return cell[z].energy; }};
    pvd.variables += {"P",   [z](Cell cell) -> double { return cell[z].pressure; }};
    pvd.variables += {"T",   [z](Cell cell) -> double { return cell[z].temperature; }};
    pvd.variables += {"b0",  [z](Cell cell) -> double { return cell[z].mass_frac[0]; }};
    pvd.variables += {"b1",  [z](Cell cell) -> double { return cell[z].mass_frac[1]; }};
    pvd.variables += {"a0",  [z](Cell cell) -> double { return cell[z].alpha(0); }};
    pvd.variables += {"a1",  [z](Cell cell) -> double { return cell[z].alpha(1); }};
    pvd.variables += {"rho0",[z](Cell cell) -> double { return cell[z].densities[0]; }};
    pvd.variables += {"rho1",[z](Cell cell) -> double { return cell[z].densities[1]; }};
    pvd.variables += {"p",  [p=data.p](Cell cell) -> double { return cell[p][0]; }};
    pvd.variables += {"n.x", [n=data.n](Cell cell) -> double { return cell[n][0].x(); }};
    pvd.variables += {"n.y", [n=data.n](Cell cell) -> double { return cell[n][0].y(); }};
    pvd.variables += {"n.z", [n=data.n](Cell cell) -> double { return cell[n][0].z(); }};

    // Initial conditions (adaptive to initial data)
    for (int k = 0; mesh.adaptive() && k < mesh.max_level() + 2; ++k) {
        init_cells(mesh, mixture, data.init);
        solver.set_flags(mesh);
        mesh.make_shuba(2);
        mesh.refine();
    }
    init_cells(mesh, mixture, data.init);

    size_t n_step = 0;
    double curr_time = 0.0;
    double next_write = 0.0;
    double max_time = 0.2;

    while (curr_time < max_time) {
        if (curr_time >= next_write) {
            std::cout << "\tStep: " << std::setw(6) << n_step << ";"
                      << "\tTime: " << std::setw(6) << std::setprecision(3) << curr_time << "\n";

            pvd.save(mesh, curr_time - 1.0);
            solver.interface_recovery(mesh);
            pvd.save(mesh, curr_time);

            auto body = solver.domain(mesh, 1);
            pvd_body.save(body, curr_time);

            next_write += max_time / 100;
        }

        // Finish exactly at max_time
        solver.set_max_dt(max_time - curr_time);

        // Integration step
        solver.update(mesh);
        solver.set_flags(mesh);
        mesh.make_shuba(2);
        mesh.refine();

        curr_time += solver.dt();
        n_step += 1;
    }
    pvd.save(mesh, max_time);

    solver.interface_recovery(mesh);
    auto body = solver.domain(mesh, 1);
    pvd_body.save(body, curr_time);

    return 0;
}
