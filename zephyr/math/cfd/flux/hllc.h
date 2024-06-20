#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>

namespace zephyr::math {

///@brief Вычисление потока методом HLLC с использованием формул для каждой величины отдельно
class HLLC : public NumFlux {
public:

    HLLC() = default;

    template<class ...Args>
    inline static std::unique_ptr<HLLC> create(Args &&... args) {
        return std::unique_ptr<HLLC>(new HLLC(std::forward<Args>(args)...));
    }

    [[nodiscard]] std::string get_name() const override { return "HLLC"; }

    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    [[nodiscard]] smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;

    [[nodiscard]] mmf::Flux
    mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const override;
};

/// @brief Оптимизированный расчёт потока по методу HLLC
class HLLC2 : public NumFlux {
public:

    HLLC2() = default;

    template<class ...Args>
    inline static std::unique_ptr<HLLC2> create(Args &&... args) {
        return std::unique_ptr<HLLC2>(new HLLC2(std::forward<Args>(args)...));
    }

    [[nodiscard]] std::string get_name() const override { return "HLLC2"; }

    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    [[nodiscard]] smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;

    [[nodiscard]] mmf::Flux
    mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const override;
};

/// @brief Central formulation of the HLLC flux
class HLLC_central : public NumFlux {
public:

    HLLC_central() = default;

    template<class ...Args>
    inline static std::unique_ptr<HLLC_central> create(Args &&... args) {
        return std::unique_ptr<HLLC_central>(new HLLC_central(std::forward<Args>(args)...));
    }

    [[nodiscard]] std::string get_name() const override { return "HLLC_central"; }

    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    [[nodiscard]] smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;

    [[nodiscard]] mmf::Flux
    mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const override;
};

/// @brief HLLC-LM flux with low Mach number correction
class HLLC_LM : public NumFlux {
public:

    HLLC_LM() = default;

    template<class ...Args>
    inline static std::unique_ptr<HLLC_LM> create(Args &&... args) {
        return std::unique_ptr<HLLC_LM>(new HLLC_LM(std::forward<Args>(args)...));
    }

    [[nodiscard]] std::string get_name() const override { return "HLLC_LM"; }

    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    [[nodiscard]] smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;

    [[nodiscard]] mmf::Flux
    mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const override;
};


}