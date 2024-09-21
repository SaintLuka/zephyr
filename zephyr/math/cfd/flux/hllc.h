#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>

namespace zephyr::math {

/// @brief Вычисление потока методом HLLC. 
/// Используется симметричная формула, можно найти в работе.
/// Nico Fleischmann, Stefan Adami, Nikolaus A.Adams. A shock-stable modification of the HLLC
/// Riemann Solver with reduced numerical dissipation. Journal of Computational Physics, 2020.
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

/// @brief HLLC-LM flux with low Mach number correction
/// Nico Fleischmann, Stefan Adami, Nikolaus A.Adams. A shock-stable modification of the HLLC
/// Riemann Solver with reduced numerical dissipation. Journal of Computational Physics, 2020.
class HLLC_LM : public NumFlux {
public:
    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<HLLC_LM>;

    /// @brief Создать умный указатель
    inline static HLLC_LM::Ptr create() {
        return std::make_shared<HLLC_LM>();
    }

    /// @brief Имя метода
    std::string get_name() const final { return "HLLC-LM"; }


    /// @brief Статическая одноматериальная версия
    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    /// @brief Одноматериальная версия
    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;
};

/// @brief HLLC-M flux (сглаживание тангенцальных компонент скорости)
/// U.S. Vevek, B. Zang, T.H. New. A carbuncle cure for the HLLC scheme using a novel
/// velocity-based sensor. Appl, Math. Mech. 2021. (есть ссылка на схему)
class HLLC_M : public NumFlux {
public:
    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<HLLC_M>;

    /// @brief Создать умный указатель
    inline static HLLC_M::Ptr create() {
        return std::make_shared<HLLC_M>();
    }

    /// @brief Имя метода
    std::string get_name() const final { return "HLLC-M"; }


    /// @brief Статическая одноматериальная версия
    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    /// @brief Одноматериальная версия
    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;
};

} // namespace zephyr::math