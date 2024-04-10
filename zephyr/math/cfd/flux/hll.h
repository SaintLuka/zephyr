#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>

namespace zephyr::math {

///@brief Вычисление потока методом HLL с использованием формул для каждой величины отдельно
class HLL : public NumFlux {
public:

    HLL() = default;

    template<class ...Args>
    inline static std::unique_ptr<HLL> create(Args &&... args) {
        return std::make_unique<HLL>(std::forward<Args>(args)...);
    }

    [[nodiscard]] std::string get_name() const override { return "HLL"; }

    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    [[nodiscard]] smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;

    [[nodiscard]] mmf::Flux
    mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const override;

    static mmf::Flux calc_mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture);

};
}