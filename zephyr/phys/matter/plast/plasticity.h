#pragma once

#include <memory>
#include <string>

namespace zephyr::phys {

/// @brief Класс для константных моделей упруго-пластичности.
/// Базовый класс для других моделей упруго-пластики.
///
/// @details
/// Постоянные модуль Юнга, модуль сдвига, коэффициент Пуассона и предел
/// текучести. Модули упругости должны быть согласованы с коэффициентом
/// Пуассона, но могут быть не согласованы с объемным модулем упругости,
/// который вычисляется отдельно в классах уравнений состояния.
class Plasticity {
public:
    double G;   ///< Модуль сдвига
    double E;   ///< Модуль Юнга
    double nu;  ///< Коэффициент Пуассона
    double Y;   ///< Предел текучести


    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<Plasticity>;
    using Ref = const std::shared_ptr<Plasticity>&;

    /// @brief Нулевые коэффициенты для модуля Юнга, модуля сдвига и предела
    /// текучести, коэффициент Пуассона равен 1/2.
    explicit Plasticity();

    /// @brief Конструктор для известных материалов
    explicit Plasticity(const std::string &name);

    /// @brief Создание указателя на базовый класс Eos
    template<typename... Args>
    static Plasticity::Ptr create(Args &&... args) {
        return std::make_shared<Plasticity>(std::forward<Args>(args)...);
    };

    /// @brief Модуль сдвига
    double shear() const;

    /// @brief Модуль Юнга
    double young() const;

    /// @brief Коэффициент Пуассона
    double poisson() const;

    /// @brief Предел текучести
    double yield() const;

protected:
    /// @brief Установить табличные параметры по названию
    void table_params(const std::string& name);
};

} // namespace zephyr::phys
