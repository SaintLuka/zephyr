#pragma once

#include <memory>

#include <zephyr/math/cfd/fluxes.h>

namespace zephyr {
namespace math {

/// @brief Расчёт потока по методу Русанова
/// @details Максимальное собственное значение считается из вектора состояния на грани
class Rusanov : public NumFlux {
public:

    Rusanov() = default;

    template<class ...Args>
    inline static std::unique_ptr<Rusanov> create(Args &&... args) {
        return std::make_unique<Rusanov>(std::forward<Args>(args)...);
    }

    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;
};

/// @brief Расчёт потока по методу русанова
/// @details Максимальное собственное значение считается из векторов состояний двух ячеек
class Rusanov2 : public NumFlux {
public:

    Rusanov2() = default;

    template<class ...Args>
    inline static std::unique_ptr<Rusanov2> create(Args &&... args) {
        return std::make_unique<Rusanov2>(std::forward<Args>(args)...);
    }

    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;
};

}
}
