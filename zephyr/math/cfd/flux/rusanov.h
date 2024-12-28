#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>

namespace zephyr::math {

/// @brief Расчёт потока методом Русанова
/// @details Предельный случай HLL при выборе |S_L| = |S_R|.
class Rusanov : public NumFlux {
public:
    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<Rusanov>;

    /// @brief Создать умный указатель
    inline static Rusanov::Ptr create() {
        return std::make_shared<Rusanov>();
    }

    /// @brief Имя метода
    std::string get_name() const final { return "Rusanov"; }


    /// @brief Статическая одноматериальная версия
    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    /// @brief Одноматериальная версия
    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;


    /// @brief Статическая многоматериальная версия
    static mmf::Flux calc_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::MixturePT &mixture);

    /// @brief Многоматериальная версия
    mmf::Flux flux(const mmf::PState &zL, const mmf::PState &zR, const phys::MixturePT &mixture) const final;
};

}
