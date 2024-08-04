#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>

namespace zephyr::math {

///@brief Вычисление потока методом HLLC
class HLLC : public NumFlux {
public:
    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<HLLC>;

    /// @brief Создать умный указатель
    inline static HLLC::Ptr create() {
        return std::make_shared<HLLC>();
    }

    /// @brief Имя метода
    std::string get_name() const final { return "HLLC"; }


    /// @brief Статическая одноматериальная версия
    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    /// @brief Одноматериальная версия
    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;


    /// @brief Статическая многоматериальная версия
    static mmf::Flux calc_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture);

    /// @brief Многоматериальная версия
    mmf::Flux flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const final;
};

/// @brief Central formulation of the HLLC flux
class HLLC_C : public NumFlux {
public:
    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<HLLC_C>;

    /// @brief Создать умный указатель
    inline static HLLC_C::Ptr create() {
        return std::make_shared<HLLC_C>();
    }

    /// @brief Имя метода
    std::string get_name() const final { return "HLLC_C"; }


    /// @brief Статическая одноматериальная версия
    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    /// @brief Одноматериальная версия
    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;

};

/// @brief HLLC-LM flux with low Mach number correction
class HLLC_LM : public NumFlux {
public:
    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<HLLC_LM>;

    /// @brief Создать умный указатель
    inline static HLLC_LM::Ptr create() {
        return std::make_shared<HLLC_LM>();
    }

    /// @brief Имя метода
    std::string get_name() const final { return "HLLC_LM"; }


    /// @brief Статическая одноматериальная версия
    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    /// @brief Одноматериальная версия
    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;
};

} // namespace zephyr::math