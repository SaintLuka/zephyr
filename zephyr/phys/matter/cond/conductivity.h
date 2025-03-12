#pragma once

#include <memory>
#include <string>

namespace zephyr::phys {

/// @brief Класс для теплопроводности. Не абстрактный, по умолчанию реализует
/// материал с нулевой теплопроводностью.
class Conductivity {
public:
    /// @brief Умный указатель на класс
    using Ptr = std::shared_ptr<Conductivity>;
    using Ref = const std::shared_ptr<Conductivity>&;

    /// @brief Создание класса с единственным параметром --
    /// коэффициентом теплопроводности
    explicit Conductivity(double kappa = 0.0);

    /// @brief Конструктор для известных материалов
    explicit Conductivity(const std::string &name);

    /// @brief Создание умного указателя на класс
    template<typename... Args>
    static Conductivity::Ptr create(Args &&... args) {
        return std::make_shared<Conductivity>(std::forward<Args>(args)...);
    };

    /// @brief Коэффициент теплопроводности
    double kappa(double temperature = 0.0) const;

protected:

    /// @brief Установить табличные параметры по названию
    void table_params(const std::string& name);


    /// @brief Коэффициент теплопроводности
    double m_kappa = 0.0;
};

} // namespace zephyr::phys
