#pragma once

#include <memory>

#include <zephyr/phys/eos/eos.h>
#include <zephyr/phys/eos/materials.h>
#include <zephyr/math/cfd/models.h>

namespace zephyr::math {

enum class Fluxes {
    CIR1,
    CIR2,
    GODUNOV,
    HLL,
    HLLC,
    HLLC2,
    RUSANOV,
    RUSANOV2
};

enum class Fluxes {
    CIR1,
    CIR2,
    GODUNOV,
    HLL,
    HLLC,
    HLLC2,
    RUSANOV,
    RUSANOV2
};

/// @brief Абстрактный класс для вычисления численного потока
class NumFlux {
public:
    using Ptr = std::unique_ptr<NumFlux>;

    /// @brief Создать нужный класс по enum
    static NumFlux::Ptr create(Fluxes flux);

    /// @brief Численный поток для классической задачи
    [[nodiscard]] virtual smf::Flux flux(const smf::PState& zL, const smf::PState& zR, const phys::Eos& eos) const {
        throw std::runtime_error("NumFlux::flux(PState...) is not implemented");
    }

    [[nodiscard]] virtual mmf::Flux mm_flux(const mmf::PState& zL, const mmf::PState& zR, const phys::Materials& mixture) const {
        throw std::runtime_error("MmNumFlux::mm_flux(PState...) is not implemented");
    }

    [[nodiscard]] virtual std::string get_name() const { return "Flux"; }

    virtual ~NumFlux() = default;
};

}