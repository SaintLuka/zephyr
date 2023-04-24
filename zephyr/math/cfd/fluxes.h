#pragma once

#include <memory>

#include <zephyr/math/cfd/models.h>
#include <zephyr/phys/eos/eos.h>

namespace zephyr { namespace math {

/// @brief Абстрактный класс для вычисления численного потока
class NumFlux {
public:
    using Ptr = std::unique_ptr<NumFlux>;

    /// @brief Численный поток для классической задачи
    virtual smf::Flux flux(const smf::PState& zL, const smf::PState& zR, const phys::Eos& eos) const {
        throw std::runtime_error("NumFlux::flux(PState...) is not implemented");
    }
};

/// @brief Вычисление потока методом Годунова
class GodunovFlux : public NumFlux {
public:

    GodunovFlux() = default;

    template <class ...Args>
    inline static std::unique_ptr<GodunovFlux> create(Args&&... args) {
        return std::make_unique<GodunovFlux>(std::forward<Args>(args)...);
    }

    /// @brief Поток как решение задачи о распаде разрыва
    static smf::Flux calc_flux(const smf::PState& zL, const smf::PState& zR, const phys::Eos& eos);

    /// @brief Поток как решение задачи о распаде разрыва
    smf::Flux flux(const smf::PState& zL, const smf::PState& zR, const phys::Eos& eos) const final;

};

}
}