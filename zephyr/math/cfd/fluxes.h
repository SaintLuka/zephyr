#pragma once

#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/eos/eos.h>

namespace zephyr { namespace math {

/// @brief Абстрактный класс для вычисления численного потока
class NumFlux {
public:
    using Ptr = std::unique_ptr<NumFlux>;

    /// @brief Численный поток для классической задачи
    virtual smf::Flux flux(const smf::PState& zL, const smf::PState& zR, const phys::Eos& eos) const {
        throw std::runtime_error("NumFlux::flux(PState...) is not implemented");
    }
};

}
}