/// @file mmf_static.cpp
/// @brief Static test for two materials

#include <iomanip>

#include <zephyr/geom/grid.h>
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
    mixture.adjust_cv({1.0_kg_m3, 1000.0_kg_m3}, 1.0_bar, 300.0);

    const PState z_air(
            1.0_kg_m3,          // density
            Vector3d::Zero(),   // velocity
            1.0_bar,            // pressure
            Fractions::Pure(0), // mass fractions
            mixture);

    const PState z_water(
            1000.0_kg_m3,       // density
            Vector3d::Zero(),   // velocity
            1.0_bar,            // pressure
            Fractions::Pure(1), // mass fractions
            mixture);

    Box box = mesh.bbox();
    Vector3d c = mesh.bbox().center();
    Vector3d n = {box.sizes().y(), box.sizes().x() / 2, 0.0};
    n.normalize();

    mesh.for_each([&](EuCell &cell) {
        // Exact section by plane
        double vol_frac0 = cell.polygon().clip_area(c, n) / cell.volume();
        double vol_frac1 = 1.0 - vol_frac0;

        if (vol_frac0 < 1.0e-12) {
            cell[z] = z_water;
            return;
        }
        if (vol_frac0 > 1.0 - 1.0e-12) {
            cell[z] = z_air;
            return;
        }

        cell[z] = PState::Mix1(mixture, {vol_frac0, vol_frac1}, {z_air, z_water});
    });
}

int main() {
    threads::off();

    // Generator of a Cartesian grid
    int nx = 5;
    int ny = 175;
    generator::Rectangle gen(0.0, 0.1*nx, 0.0, 0.1*ny);
    gen.set_sizes(nx, ny);
    gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
                        .bottom=Boundary::WALL, .top=Boundary::WALL});
    gen.set_adaptive(false);

    Grid grid = gen.make();

    // Distort the grid (add random offsets to nodes)
    grid.transform(
        [](const Vector3d& v) -> Vector3d {
            Vector3d res = v;
            res.x() += 0.002 * (2.0 * rand() / double(RAND_MAX) - 1.0);
            res.y() += 0.002 * (2.0 * rand() / double(RAND_MAX) - 1.0);
            return res;
        });
    grid.make_amr();

    // Create mesh
    EuMesh mesh(std::move(grid));

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
    solver.set_splitting(DirSplit::NONE);

    // Add data fields, choose main data layer
    auto data = solver.add_types(mesh);
    auto z = data.init;

    // Files for output
    PvdFile pvd("mesh", "output");

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
    //pvd.variables += {"e0",[z,mixture](EuCell cell) -> double { return cell[z].true_energy(mixture, 0); }};
    //pvd.variables += {"e1",[z,mixture](EuCell cell) -> double { return cell[z].true_energy(mixture, 1); }};
    //pvd.variables += {"n.x", [n=data.n](EuCell cell) -> double { return cell[n][0].x(); }};
    //pvd.variables += {"n.y", [n=data.n](EuCell cell) -> double { return cell[n][0].y(); }};

    // Initial conditions
    init_cells(mesh, mixture, data.init);

    size_t n_step = 0;
    double curr_time = 0.0;

    pvd.save(mesh, 0.0);
    while (n_step < 100'000) {
        std::cout << "\tStep: " << std::setw(6) << n_step << ";\n";

        if (n_step % 1000 == 0) {
            pvd.save(mesh, n_step);
        }

        // Integration step
        solver.update(mesh);
        curr_time += solver.dt();
        n_step += 1;
    }

    return 0;
}
