#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>

namespace zephyr::math {

/// @brief Вычисление потока методом Годунова
class Godunov : public NumFlux{
public:

    Godunov() = default;

    template<class ...Args>
    inline static std::unique_ptr<Godunov> create(Args &&... args) {
        return std::unique_ptr<Godunov>(new Godunov(std::forward<Args>(args)...));
    }

    [[nodiscard]] std::string get_name() const override { return "Godunov"; }

    /// @brief Поток как решение задачи о распаде разрыва
    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    /// @brief Поток как решение задачи о распаде разрыва
    [[nodiscard]] smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;

    /// @brief Поток как решение задачи о распаде разрыва
    static mmf::Flux calc_mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture);

    /// @brief Поток как решение задачи о распаде разрыва
    [[nodiscard]] mmf::Flux
    mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const final;
};

}