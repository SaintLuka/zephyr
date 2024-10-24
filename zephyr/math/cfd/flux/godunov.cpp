#include <zephyr/math/solver/riemann.h>
#include <zephyr/math/cfd/flux/godunov.h>
#include <zephyr/math/cfd/flux/rusanov.h>
#include <zephyr/math/cfd/flux/hll.h>
#include <iostream>
#include <boost/format.hpp>

namespace zephyr::math {

using phys::Eos;

smf::Flux Godunov::flux(const smf::PState &zL, const smf::PState &zR, const Eos &eos) const {
    return calc_flux(zL, zR, eos);
}

smf::Flux Godunov::calc_flux(const smf::PState &zL, const smf::PState &zR, const Eos &eos) {
    auto sol = RiemannSolver::solve(zL, zR,
                                    eos.stiffened_gas(zL.density, zL.pressure),
                                    eos.stiffened_gas(zR.density, zR.pressure));

    Vector3d v = {
            sol.Uf,
            sol.U > 0.0 ? zL.velocity.y() : zR.velocity.y(),
            sol.U > 0.0 ? zL.velocity.z() : zR.velocity.z()
    };

    smf::PState z(sol.rho, v, sol.Pf, eos.energy_rP(sol.rho, sol.Pf));

    smf::Flux flux(z);

    return flux;
}

mmf::Flux Godunov::flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const {
    return calc_flux(zL, zR, mixture);
}

mmf::Flux Godunov::calc_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) {
    bool zL_bad = std::isinf(zL.density) || std::isnan(zL.density) ||
                  std::isinf(zL.velocity.x()) || std::isnan(zL.velocity.x()) ||
                  std::isinf(zL.velocity.y()) || std::isnan(zL.velocity.y()) ||
                  std::isinf(zL.velocity.z()) || std::isnan(zL.velocity.z()) ||
                  std::isinf(zL.pressure) || std::isnan(zL.pressure) ||
                  zL.mass_frac.empty();
    if (zL_bad) {
        std::cerr << "Bad left state in input of Godunov flux: " << zL << "\n";
        exit(1);
        throw std::runtime_error("bad left state");
    }
    bool zR_bad = std::isinf(zR.density) || std::isnan(zR.density) ||
                  std::isinf(zR.velocity.x()) || std::isnan(zR.velocity.x()) ||
                  std::isinf(zR.velocity.y()) || std::isnan(zR.velocity.y()) ||
                  std::isinf(zR.velocity.z()) || std::isnan(zR.velocity.z()) ||
                  std::isinf(zR.pressure) || std::isnan(zR.pressure) ||
                  zR.mass_frac.empty();
    if (zR_bad) {
        std::cerr << "Bad right state in input of Godunov flux: " << zR << "\n";
        exit(1);
        throw std::runtime_error("bad right state");
    }

    auto sol = RiemannSolver::solve(zL.to_smf(), zR.to_smf(),
                                    mixture.stiffened_gas(zL.density, zL.pressure, zL.mass_frac, {.T0 = zL.temperature}),
                                    mixture.stiffened_gas(zR.density, zR.pressure, zR.mass_frac, {.T0 = zR.temperature}));
    if (sol.U * sol.Uf < 0 && abs(sol.U - sol.Uf) > 0.1) {
        std::cerr << "Different sign of U and Uf\n" << "U: " << sol.U << ", Uf: " << sol.Uf << '\n';
    }

    mmf::Flux flux;
    if (sol.U >= 0) {
        Vector3d v(sol.Uf, zL.velocity.y(), zL.velocity.z());
        auto [energy, temperature] = mixture.find_eT(sol.rho, sol.Pf, zL.mass_frac, {.T0=zL.temperature});
        mmf::PState z(sol.rho, v, sol.Pf, energy, temperature, zL.mass_frac, mmf::Fractions::NaN());

        flux = mmf::Flux(z);
    } else {
        Vector3d v(sol.Uf, zR.velocity.y(), zR.velocity.z());
        auto [energy, temperature] = mixture.find_eT(sol.rho, sol.Pf, zR.mass_frac, {.T0=zR.temperature});
        mmf::PState z(sol.rho, v, sol.Pf, energy, temperature, zR.mass_frac, mmf::Fractions::NaN());

        flux = mmf::Flux(z);
    }

    if (flux.is_bad()) {
        std::cerr << "calc Godunov flux\n";
//        std::cerr << "zL: " << zL << "\nzR: " << zR << "\n";
        std::cerr << "Flux: " << flux << "\n";
        std::cerr << "Code for debug:\n";
        std::cerr << boost::format("PState zL(%1%, {%2%, %3%, %4%}, %5%, %6%, %7%, %8%);\n") %
                     zL.density % zL.velocity.x() % zL.velocity.y() % zL.velocity.z() % zL.pressure %
                     zL.energy % zL.temperature % zL.mass_frac;
        std::cerr << boost::format("PState zR(%1%, {%2%, %3%, %4%}, %5%, %6%, %7%, %8%);\n") %
                     zR.density % zR.velocity.x() % zR.velocity.y() % zR.velocity.z() % zR.pressure %
                     zR.energy % zR.temperature % zR.mass_frac;
//        std::cerr << "Godunov flux RIP\n";
//        flux = Rusanov2::calc_mm_flux(zL, zR, mixture);
        exit(1);
//        throw std::runtime_error("Bad Godunov flux");
    }

    return flux;
}

} // namespace zephyr::math