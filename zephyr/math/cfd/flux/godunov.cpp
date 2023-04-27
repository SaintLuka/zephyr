#include <zephyr/math/solver/riemann.h>
#include <zephyr/math/cfd/flux/godunov.h>

namespace zephyr { namespace math {

using smf::PState;
using smf::Flux;
using phys::Eos;

Flux Godunov::calc_flux(const PState& zL, const PState& zR, const Eos& eos) {
    auto sol = RiemannSolver::solve(zL, zR, eos);

    Vector3d v = {
            sol.U,
            sol.U > 0.0 ? zL.velocity.y() : zR.velocity.y(),
            sol.U > 0.0 ? zL.velocity.z() : zR.velocity.z()
    };

    PState z(sol.rho, v, sol.P, eos.energy_rp(sol.rho, sol.P));

    Flux flux(z);

    return flux;
}

Flux Godunov::flux(const PState& zL, const PState& zR, const Eos& eos) const {
    return calc_flux(zL, zR, eos);
}

} // namespace math
} // namespace zephyr