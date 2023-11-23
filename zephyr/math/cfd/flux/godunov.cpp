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

mmf::Flux Godunov::calc_mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) {
    auto sol = RiemannSolver::solve(zL.to_smf(), zR.to_smf(),
                                    mixture.stiffened_gas(zL.density, zL.pressure, zL.mass_frac),
                                    mixture.stiffened_gas(zR.density, zR.pressure, zR.mass_frac));

    if (sol.U >= 0) {
        Vector3d v(sol.U, zL.velocity.y(), zL.velocity.z());
        double energy = mixture.energy_rp(sol.rho, sol.P, zL.mass_frac, {.T0=zL.temperature});
        double temperature = mixture.temperature_rp(sol.rho, sol.P, zL.mass_frac, {.T0=zL.temperature});
        mmf::PState z(sol.rho, v, sol.P, energy, temperature, zL.mass_frac);

        return mmf::Flux(z);
    } else {
        Vector3d v(sol.U, zR.velocity.y(), zR.velocity.z());
        double energy = mixture.energy_rp(sol.rho, sol.P, zR.mass_frac, {.T0=zR.temperature});
        double temperature = mixture.temperature_rp(sol.rho, sol.P, zR.mass_frac, {.T0=zR.temperature});
        mmf::PState z(sol.rho, v, sol.P, energy, temperature, zR.mass_frac);

        return mmf::Flux(z);
    }
}

mmf::Flux Godunov::mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const {
    return calc_mm_flux(zL, zR, mixture);
}

} // namespace zephyr