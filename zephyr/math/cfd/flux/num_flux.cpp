#include <zephyr/math/cfd/fluxes.h>

namespace zephyr { namespace math {

NumFlux::Ptr NumFlux::create(Fluxes flux) {
    switch (flux) {
        case Fluxes::CIR1:
            return CIR1::create();
        case Fluxes::CIR2:
            return CIR2::create();
        case Fluxes::GODUNOV:
            return Godunov::create();
        case Fluxes::HLLC:
            return HLLC::create();
        case Fluxes::HLLC2:
            return HLLC2::create();
        case Fluxes::RUSANOV:
            return Rusanov::create();
        case Fluxes::RUSANOV2:
            return Rusanov2::create();
    }
}

}
}