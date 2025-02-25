#pragma once

#include <zephyr/math/cfd/flux/num_flux.h>

namespace zephyr::math {

namespace smf {
/// @brief Двухволновая конфигурация
struct WaveConfig2 {
    double S_L, S_R;  ///< Оценки минимального и максимального СЗ
    smf::QState Qs;   ///< Консервативный вектор в возмущенном регионе
    smf::Flux   Fs;   ///< Вектор потока в возмущенном регионе
};
}

namespace mmf {
/// @brief Двухволновая конфигурация
struct WaveConfig2 {
    double S_L, S_R;  ///< Оценки минимального и максимального СЗ
    mmf::QState Qs;   ///< Консервативный вектор в возмущенном регионе
    mmf::Flux   Fs;   ///< Вектор потока в возмущенном регионе
};
}

///@brief Вычисление потока методом HLL с использованием формул для каждой величины отдельно
class HLL : public NumFlux {
public:
    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<HLL>;

    /// @brief Создать умный указатель
    inline static HLL::Ptr create() {
        return std::make_shared<HLL>();
    }

    /// @brief Имя метода
    std::string get_name() const final { return "HLL"; }

    // ========================================================================
    //                    Одноматериальные версии функций
    // ========================================================================

    /// @brief Статическая одноматериальная версия
    static smf::Flux calc_flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos);

    /// @brief Одноматериальная версия
    smf::Flux flux(const smf::PState &zL, const smf::PState &zR, const phys::Eos &eos) const final;

    /// @brief Волновая конфигурация
    static smf::WaveConfig2 wave_config(const phys::Eos &eos,
                                        const smf::PState &zL, const smf::PState &zR);

    /// @brief Волновая конфигурация
    /// Потоки fL/fR не обязательно согласованы с qL/qR
    static smf::WaveConfig2 wave_config(const phys::Eos &eos,
                                        const smf::QState &qL, const smf::Flux &fL,
                                        const smf::QState &qR, const smf::Flux &fR);

    // ========================================================================
    //                    Многоматериальные версии функций
    // ========================================================================

    /// @brief Статическая многоматериальная версия
    static mmf::Flux calc_flux(const mmf::PState &zL, const mmf::PState &zR, const phys::MixturePT &mix);

    /// @brief Многоматериальная версия
    mmf::Flux flux(const mmf::PState &zL, const mmf::PState &zR, const phys::MixturePT &mix) const final;

    /// @brief Волновая конфигурация
    static mmf::WaveConfig2 wave_config(const phys::MixturePT &mix,
                                        const mmf::PState &zL, const mmf::PState &zR);

    /// @brief Волновая конфигурация
    /// Потоки fL/fR не обязательно согласованы с qL/qR
    static mmf::WaveConfig2 wave_config(const phys::MixturePT &mix,
                                        const mmf::QState &qL, const mmf::Flux &fL,
                                        const mmf::QState &qR, const mmf::Flux &fR);

};

} // namespace zephyr::math