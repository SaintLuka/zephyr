#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>
#include <zephyr/math/cfd/flux/cir.h>
#include <zephyr/math/cfd/flux/godunov.h>
#include <zephyr/math/cfd/flux/hllc.h>
#include <zephyr/math/cfd/flux/rusanov.h>


enum class Fluxes{
    CIR1,
    CIR2,
    GODUNOV,
    HLLC,
    HLLC2,
    RUSANOV,
    RUSANOV2
};

zephyr::math::NumFlux::Ptr flux_from_enum(Fluxes flux){
    switch (flux) {
        case Fluxes::CIR1:
            return zephyr::math::CIR1::create();
        case Fluxes::CIR2:
            return zephyr::math::CIR2::create();
        case Fluxes::GODUNOV:
            return zephyr::math::Godunov::create();
        case Fluxes::HLLC:
            return zephyr::math::HLLC::create();
        case Fluxes::HLLC2:
            return zephyr::math::HLLC2::create();
        case Fluxes::RUSANOV:
            return zephyr::math::Rusanov::create();
        case Fluxes::RUSANOV2:
            return zephyr::math::Rusanov2::create();
    }
}