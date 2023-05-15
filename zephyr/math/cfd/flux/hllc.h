#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>

namespace zephyr {
namespace math {


///@brief Вычисление потока методом HLLC с использованием формул для каждой величины отдельно
class HLLC : public NumFlux {
public:

    HLLC() = default;

    template<class ...Args>
    inline static std::unique_ptr<HLLC> create(Args &&... args) {
        return std::make_unique<HLLC>(std::forward<Args>(args)...);
    }

    std::string get_name() const override {return "HLLC"; }

    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;
};

/// @brief Оптимизированный расчёт потока по методу HLLC
class HLLC2 : public NumFlux {
public:

    HLLC2() = default;

    template<class ...Args>
    inline static std::unique_ptr<HLLC2> create(Args &&... args) {
        return std::make_unique<HLLC2>(std::forward<Args>(args)...);
    }

    std::string get_name() const override {return "HLLC2"; }

    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;
};

}
}