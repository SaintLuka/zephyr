#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>

namespace zephyr {
namespace math {


///@brief Вычисление потока методом HLL с использованием формул для каждой величины отдельно
class HLL : public NumFlux {
public:

    HLL() = default;

    template<class ...Args>
    inline static std::unique_ptr<HLL> create(Args &&... args) {
        return std::make_unique<HLL>(std::forward<Args>(args)...);
    }

    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;
};
}
}