#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>

namespace zephyr { namespace math {

/// @brief Вычисление потока методом Годунова
class Godunov : public NumFlux {
public:

    Godunov() = default;

    template <class ...Args>
    inline static std::unique_ptr<Godunov> create(Args&&... args) {
        return std::make_unique<Godunov>(std::forward<Args>(args)...);
    }

    std::string get_name() const override {return "Godunov"; }

    /// @brief Поток как решение задачи о распаде разрыва
    static smf::Flux calc_flux(const smf::PState& zL, const smf::PState& zR, const phys::Eos& eos);

    /// @brief Поток как решение задачи о распаде разрыва
    smf::Flux flux(const smf::PState& zL, const smf::PState& zR, const phys::Eos& eos) const final;

};

}
}