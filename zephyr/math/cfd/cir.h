#pragma once

#include <memory>

#include <zephyr/math/cfd/fluxes.h>

namespace zephyr { namespace math {

/// @brief Вычисление потока по схеме Куранта-Изаксона-Риса,
/// классическая схема
class CIR1 : public NumFlux {
public:

    CIR1() = default;

    template <class ...Args>
    inline static std::unique_ptr<CIR1> create(Args&&... args) {
        return std::make_unique<CIR1>(std::forward<Args>(args)...);
    }

    static smf::Flux calc_flux(const smf::PState& zL, const smf::PState& zR, const phys::Eos& eos);

    smf::Flux flux(const smf::PState& zL, const smf::PState& zR, const phys::Eos& eos) const final;

};

/// @brief Вычисление потока по схеме Куранта-Изаксона-Риса,
/// схема без использования производных от уравнений состояния
class CIR2 : public NumFlux {
public:

    CIR2() = default;

    template <class ...Args>
    inline static std::unique_ptr<CIR2> create(Args&&... args) {
        return std::make_unique<CIR2>(std::forward<Args>(args)...);
    }

    static smf::Flux calc_flux(const smf::PState& zL, const smf::PState& zR, const phys::Eos& eos);

    smf::Flux flux(const smf::PState& zL, const smf::PState& zR, const phys::Eos& eos) const final;

};

}
}
