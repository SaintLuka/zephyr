#include <zephyr/math/cfd/fluxes.h>
#include "num_flux.h"


namespace zephyr::math {

NumFlux::Ptr NumFlux::create(Fluxes flux) {
    switch (flux) {
        case Fluxes::CIR1:
            return CIR1::create();
        case Fluxes::CIR2:
            return CIR2::create();
        case Fluxes::GODUNOV:
            return Godunov::create();
        case Fluxes::HLL:
            return HLL::create();
        case Fluxes::HLLC:
            return HLLC::create();
        case Fluxes::HLLC2:
            return HLLC2::create();
        case Fluxes::RUSANOV:
            return Rusanov::create();
        case Fluxes::RUSANOV2:
            return Rusanov2::create();
        case Fluxes::HLLC_central:
            return HLLC_central::create();
        case Fluxes::HLLC_LM:
            return HLLC_LM::create();
        default:
            throw std::runtime_error("Unknown numerical flux");
    }
}
}