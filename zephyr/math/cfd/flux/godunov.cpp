#include <zephyr/math/solver/riemann.h>
#include <zephyr/math/cfd/flux/godunov.h>

namespace zephyr::math {

using phys::Eos;

smf::Flux Godunov::calc_flux(const smf::PState &zL, const smf::PState &zR, const Eos &eos) {
    auto sol = RiemannSolver::solve(zL, zR,
                                    eos.stiffened_gas(zL.density, zL.pressure),
                                    eos.stiffened_gas(zR.density, zR.pressure));

    Vector3d v = {
            sol.U,
            sol.U > 0.0 ? zL.velocity.y() : zR.velocity.y(),
            sol.U > 0.0 ? zL.velocity.z() : zR.velocity.z()
    };

    smf::PState z(sol.rho, v, sol.P, eos.energy_rp(sol.rho, sol.P));

    smf::Flux flux(z);

    return flux;
}

smf::Flux Godunov::flux(const smf::PState &zL, const smf::PState &zR, const Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

mmf::Flux Godunov::calc_mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Eos &eos) {
    auto sol = RiemannSolver::solve(zL.to_smf(), zR.to_smf(),
                                    eos.stiffened_gas(zL.density, zL.pressure),
                                    eos.stiffened_gas(zR.density, zR.pressure));

    Vector3d v = {
            sol.U,
            sol.U > 0.0 ? zL.velocity.y() : zR.velocity.y(),
            sol.U > 0.0 ? zL.velocity.z() : zR.velocity.z()
    };
    double energy = eos.energy_rp(sol.rho, sol.P);
    double temperature = eos.temperature_rp(sol.rho, sol.P);

    if (sol.U > 0) {
        mmf::PState z(sol.rho, v, sol.P, energy, temperature, zL.mass_frac);
        return mmf::Flux(z);
    } else {
        mmf::PState z(sol.rho, v, sol.P, energy, temperature, zR.mass_frac);
        return mmf::Flux(z);
    }
}

mmf::Flux Godunov::mm_flux(const mmf::PState &zL, const mmf::PState &zR, const Eos &eos) const {
    return calc_mm_flux(zL, zR, eos);
}

} // namespace zephyr