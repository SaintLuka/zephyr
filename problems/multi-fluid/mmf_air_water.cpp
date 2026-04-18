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

#include <zephyr/geom/sections.h>
#include <zephyr/math/calc/roots.h>
#include <zephyr/utils/numpy.h>
#include <zephyr/utils/pyplot.h>

using namespace zephyr;
using namespace zephyr::phys;
using namespace zephyr::math;
using namespace zephyr::mesh;
using namespace zephyr::math::mmf;

using zephyr::io::PvdFile;
using zephyr::utils::mpi;
using zephyr::utils::threads;
using zephyr::utils::Stopwatch;


double temperature_dichotomy(MixturePT& mix, double P, double e, const Fractions& a) {

    auto func = [&](double T) -> double {
        double res = 0.0;
        for (int i = 0; i < 2; ++i) {
            if (a[i] > 0.0) {
                res += a[i] * (mix[i].energy_PT(P, T) - e) / mix[i].volume_PT(P, T);
            }
        }
        return res;
    };

    //std::vector Ts = np::linspace(231.0, 360.0, 1000);
    //std::vector Fs = np::zeros_like(Ts);
    //for (int i = 0; i < Ts.size(); ++i) {
    //    Fs[i] = func(Ts[i]); //(mix[1].energy_PT(P, Ts[i]) - e) / mix[1].volume_PT(P, Ts[i]);
    //}
    //utils::pyplot plt;
    //plt.plot(Ts, Fs);
    //plt.show();

    return dichotomy(func, 231.0, 360.0);
}

double temperature_newton(MixturePT& mix, double P, const ScalarSet& T, const Fractions& a) {
    double F = 0.0;
    for (int i = 0; i < mix.size(); ++i) {
        if (a.has(i)) {
            F += a[i] * mix[i].energy_PT(P, T[i]) / mix[i].volume_PT(P, T[i]);
        }
    }

    auto func = [mix, P, a, F](double T) -> double {
        double res = 0.0;
        for (int i = 0; i < mix.size(); ++i) {
            if (a.has(i)) {
                res += a[i] * mix[i].energy_PT(P, T) / mix[i].volume_PT(P, T);
            }
        }
        return res - F;
    };

    auto deriv = [mix, P, a](double T) -> double {
        double res = 0.0;
        for (int i = 0; i < mix.size(); ++i) {
            if (a.has(i)) {
                auto v = mix[i].volume_PT(P, T, {.deriv=true});
                auto e = mix[i].energy_PT(P, T, {.deriv=true});
                res += a[i] * (e.dT / v - e * v.dT / v / v);
            }
        }
        return res;
    };

    double T_min = std::min(T[0], T[1]);
    double T_max = std::max(T[0], T[1]);

    double Tref = a[0] * T[0] + a[1] * T[1];

    double res = newton_method(func, deriv, Tref, {.max_iters = 1, .x_min = T_min, .x_max = T_max});

    if (std::isnan(res) || std::isinf(res)) {
        std::vector Ts = np::linspace(T_min, T_max, 1000);
        std::vector Fs = np::zeros_like(Ts);
        for (int i = 0; i < Ts.size(); ++i) {
            Fs[i] = func(Ts[i]); //(mix[1].energy_PT(P, Ts[i]) - e) / mix[1].volume_PT(P, Ts[i]);
        }
        utils::pyplot plt;
        plt.plot(Ts, Fs);
        //plt.show();
        
        double e0 = mix[0].energy_PT(P, T[0]);
        double e1 = mix[1].energy_PT(P, T[1]);
        
        double e0ref = mix[0].energy_PT(P, Tref);
        double e1ref = mix[1].energy_PT(P, Tref);
        
        double v0 = mix[0].volume_PT(P, T[0]);
        double v1 = mix[1].volume_PT(P, T[1]);
        
        double v0ref = mix[0].volume_PT(P, Tref);
        double v1ref = mix[1].volume_PT(P, Tref);
        
        

        double ololo = func(Tref);

        newton_method<true>(func, deriv, Tref, {.max_iters = 1, .x_min = T_min, .x_max = T_max});
    }
    return res;
}

double temperature_min_energy(MixturePT& mix, double P, const ScalarSet& T, const Fractions& a) {
    double T_min = std::min(T[0], T[1]);
    double T_max = std::max(T[0], T[1]);

    if (a[0] == 0.0) { return T[1]; }
    if (a[1] == 0.0) { return T[0]; }

    std::vector Ts = np::linspace(T_min, T_max, 100);

    auto func = [mix, P, a](double T) -> double {
        double v0 = mix[0].volume_PT(P, T);
        double v1 = mix[1].volume_PT(P, T);
        double v = 1.0 / (a[0] / v0 + a[1] / v1);
        double b0 = a[0] * v / v0;
        double b1 = a[1] * v / v1;
        double e = b0 * mix[0].energy_PT(P, T) + b1 * mix[1].energy_PT(P, T);
        return e;
    };

    double e_min = 1.0e300;
    double Tres = -1.0;
    for (auto t: Ts) {
        double e = func(t);
        if (e < e_min) {
            e_min = e;
            Tres = t;
        }
    }

    /*
    std::vector Fs = np::zeros_like(Ts);
    for (int i = 0; i < Ts.size(); ++i) {
        Fs[i] = func(Ts[i]);
    }
    utils::pyplot plt;
    plt.plot(Ts, Fs);
    plt.show();*/

    return Tres;
}

void init_cells(EuMesh& mesh, MixturePT& mixture, Storable<PState> z) {
    //mixture.adjust_cv({1.0_kg_m3, 1000.0_kg_m3}, 1.0_bar, 300.0);

    const PState za(
            1.0_kg_m3,          // density
            Vector3d::Zero(),   // velocity
            1.0_bar,            // pressure
            Fractions::Pure(0), // mass fractions
            mixture);

    const PState zw(
            1000.0_kg_m3,       // density
            Vector3d::Zero(),   // velocity
            1.0_bar,            // pressure
            Fractions::Pure(1), // mass fractions
            mixture);

    for (int i = 0; i < mesh.n_cells(); ++i) {
        EuCell cell = mesh[i];

        const double eps = 1.0e-9;
        double x1 = cell.face(Side2D::L).x();
        double x2 = cell.face(Side2D::R).x();

        /*
        // Слева вода
        if (x2 < eps) {
            cell[z] = zw;
            continue;
        }
        // Справа газ
        if (x1 > eps) {
            cell[z] = za;
            continue;
        }
        */

        // Mixed cell (air and water)
        const Vector3d normal = {0.9, -0.1, 0.0};
        double vol_frac_air = quad_volume_fraction(cell.center().dot(normal), normal, cell.hx(), cell.hy());
        
        PState Za = za;
        PState Zw = zw;

        //double T = vol_frac_air * Za.temperature + (1.0 - vol_frac_air) * Zw.temperature;
        //double e = vol_frac_air * za.energy + (1.0 - vol_frac_air) * zw.energy;
        //double T = temperature_dichotomy(mixture, za.P(), e, {vol_frac_air, 1.0 - vol_frac_air});
        double T = temperature_min_energy(mixture, za.P(), {za.T(), zw.T()}, {vol_frac_air, 1.0 - vol_frac_air});
        Za.density = 1.0 / mixture[0].volume_PT(Za.pressure, T);
        Zw.density = 1.0 / mixture[1].volume_PT(Zw.pressure, T);

        // rho = sum a_i rho_i
        double    density   = vol_frac_air * Za.density + (1.0 - vol_frac_air) * Zw.density;
        Fractions mass_frac = {vol_frac_air * Za.density / density, (1.0 - vol_frac_air) * Zw.density / density};
        mass_frac.normalize();
        Vector3d  velocity  = mass_frac[0] * Za.velocity + mass_frac[1] * Zw.velocity;
        double    pressure  = vol_frac_air * Za.pressure + (1.0 - vol_frac_air) * Zw.pressure;

        PState mix(density, velocity, pressure, mass_frac, mixture);

        cell[z] = mix;
    }
}

int main(int argc, char** argv) {
    mpi::handler handler(argc, argv);
    threads::init(argc, argv);
    threads::info();
    threads::off();

    // Generator of a Cartesian grid
    generator::Rectangle gen(-0.1_cm, 0.1_cm, 0.0_cm, 0.1_cm);
    gen.set_nx(40);
    gen.set_boundaries({.left=Boundary::ZOE, .right=Boundary::ZOE,
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
    PvdFile pvd("mesh", "output");
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
    //pvd.variables += {"e0",[z,mixture](EuCell cell) -> double { return cell[z].true_energy(mixture, 0); }};
    //pvd.variables += {"e1",[z,mixture](EuCell cell) -> double { return cell[z].true_energy(mixture, 1); }};
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
    double max_time = 4.5_us;

    Stopwatch elapsed(true);
    while (curr_time < max_time) {
        if (curr_time >= next_write) {
            //std::cout << "\tStep: " << std::setw(6) << n_step << ";"
            //          << "\tTime: " << std::setw(6) << std::setprecision(3) << 1.0e6 * curr_time << " us\n";
            pvd.save(mesh, curr_time);

            solver.interface_recovery(mesh);
            auto bubble = solver.domain(mesh, 0);
            pvd_bubble.save(bubble, curr_time);

            next_write += 0.0;//max_time / 50;
        }

        // Finish exactly at max_time
        // solver.set_max_dt(max_time - curr_time);

        // Integration step
        solver.update(mesh);
        //for (auto cell: mesh) { cell.set_flag(1); }
        solver.set_flags(mesh);
        mesh.refine();

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
