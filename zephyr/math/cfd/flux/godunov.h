#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>

namespace zephyr::math {

/// @brief Вычисление потока методом Годунова
class Godunov : public NumFlux{
public:
    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<Godunov>;

    /// @brief Создать умный указатель
    inline static Godunov::Ptr create() {
        return std::make_shared<Godunov>();
    }

    /// @brief Имя метода
    std::string get_name() const final { return "Godunov"; }


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