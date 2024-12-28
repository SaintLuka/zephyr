#pragma once

#include <memory>

#include <zephyr/phys/matter/eos/eos.h>
#include <zephyr/phys/matter/mixture_pt.h>
#include <zephyr/math/cfd/models.h>

namespace zephyr::math {

enum class Fluxes {
    GODUNOV,
    RUSANOV,
    HLL,
    HLLC,
    HLLC_M,
    HLLC_LM,
    CIR1,
    CIR2,
    CRP    ///< Композитная задача Римана, для многоматериальных
};

/// @brief Абстрактный класс для вычисления численного потока
class NumFlux {
public:
    using Ptr = std::shared_ptr<NumFlux>;

    /// @brief Создать нужный класс по enum
    static NumFlux::Ptr create(Fluxes flux);

    /// @brief Виртуальный деструктор
    virtual ~NumFlux() = default;

    /// @brief Название расчетного метода
    virtual std::string get_name() const { return "Flux"; }

    /// @brief Численный поток для одноматериальной задачи
    virtual smf::Flux flux(const smf::PState& zL, const smf::PState& zR, const phys::Eos& eos) const {
        throw std::runtime_error("NumFlux::flux(smf::PState...) is not implemented");
    }

    /// @brief Численный поток для многоматериальной задачи
    virtual mmf::Flux flux(const mmf::PState& zL, const mmf::PState& zR, const phys::MixturePT& mixture) const {
        throw std::runtime_error("NumFlux::flux(mmf::PState...) is not implemented");
    }

};

} // namespace zephyr::math