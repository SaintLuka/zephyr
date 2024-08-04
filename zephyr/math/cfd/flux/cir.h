#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>

namespace zephyr::math {

/// @brief Вычисление потока по схеме Куранта-Изаксона-Риса,
/// классическая схема
class CIR1 : public NumFlux {
public:
    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<CIR1>;

    /// @brief Создать умный указатель
    inline static CIR1::Ptr create() {
        return std::make_shared<CIR1>();
    }

    /// @brief Имя метода
    std::string get_name() const final { return "CIR1"; }


    /// @brief Статическая одноматериальная версия
    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    /// @brief Одноматериальная версия
    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;

};

/// @brief Вычисление потока по схеме Куранта-Изаксона-Риса,
/// схема без использования производных от уравнений состояния
class CIR2 : public NumFlux {
public:
    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<CIR2>;

    /// @brief Создать умный указатель
    inline static CIR2::Ptr create() {
        return std::make_shared<CIR2>();
    }

    /// @brief Имя метода
    std::string get_name() const final { return "CIR2"; }


    /// @brief Статическая одноматериальная версия
    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    /// @brief Одноматериальная версия
    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;

};

} // namespace zephyr::math
