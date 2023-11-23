#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>

namespace zephyr::math {

/// @brief Расчёт потока по методу Русанова
/// @details Максимальное собственное значение считается из вектора состояния на грани
class Rusanov : public NumFlux {
public:

    Rusanov() = default;

    template<class ...Args>
    inline static std::unique_ptr<Rusanov> create(Args &&... args) {
        return std::unique_ptr<Rusanov>(new Rusanov(std::forward<Args>(args)...));
    }

    [[nodiscard]] std::string get_name() const override {return "Rusanov"; }

    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    [[nodiscard]] smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;
};

/// @brief Расчёт потока по методу русанова
/// @details Максимальное собственное значение считается из векторов состояний двух ячеек
class Rusanov2 : public NumFlux {
public:

    Rusanov2() = default;

    template<class ...Args>
    inline static std::unique_ptr<Rusanov2> create(Args &&... args) {
        return std::unique_ptr<Rusanov2>(new Rusanov2(std::forward<Args>(args)...));
    }

    [[nodiscard]] std::string get_name() const override {return "Rusanov2"; }

    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    [[nodiscard]] smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;

    mmf::Flux mm_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::Materials &mixture) const override;
};

}
