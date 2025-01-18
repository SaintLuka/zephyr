#pragma once

#include <memory>
#include <string>

namespace zephyr::phys {

/// @class Класс для константных моделей вязкости, а так же
/// базовый класс для других моделей вязкости.
class Viscosity {
public:
    double nu;   ///< Кинематическая вязкость
    double eta;  ///< Сдвиговая/динамическая вязкость
    double zeta; ///< Объемная/вторая вязкость


    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<Viscosity>;
    using Ref = const std::shared_ptr<Viscosity>&;

    /// @brief Вязкость отсутствует, все коэффициенты равны нулю.
    explicit Viscosity();

    /// @brief Конструктор с заданием кинематической вязкости,
    /// для несжимаемых жидкостей
    explicit Viscosity(double nu);

    /// @brief Конструктор с заданием сдвиговой и объемной
    /// вязкостей, для сжимаемых жидкостей
    explicit Viscosity(double eta, double zeta);

    /// @brief Конструктор для известных материалов
    explicit Viscosity(const std::string &name);

    /// @brief Создание указателя на базовый класс Eos
    template<typename... Args>
    static Viscosity::Ptr create(Args &&... args) {
        return std::make_shared<Viscosity>(std::forward<Args>(args)...);
    };

    /// @brief Кинематическая вязкость
    double kinematic_visc(double temperature = 0.0) const;

    /// @brief Свдиговая/динамическая вязкость
    double shear_visc(double temperature = 0.0) const;

    /// @brief Объемная/вторая вязкость
    double volume_visc(double temperature = 0.0) const;

protected:
    /// @brief Установить табличные параметры по названию
    void table_params(const std::string& name);
};

} // namespace zephyr::phys
