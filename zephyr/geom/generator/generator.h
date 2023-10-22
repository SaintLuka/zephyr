#pragma once

#include <string>
#include <memory>

namespace zephyr::geom {

class Box;  ///< Ограничивающий объем
class Grid; ///< Сетка общего вида, создается генератором

/// @class Базовый класс для сеточных генераторов. Содержит набор
/// виртуальных функций для создания Storage с сеткой.
class Generator {
public:
    using Ptr = std::shared_ptr<Generator>;

    /// @brief Базовый конструктор только с именем типа сетки
    explicit Generator(const std::string &name);

    /// @brief Тип сеточного генератора
    const std::string &name() const;

    /// @brief Количество ячеек сетки
    virtual int size() const = 0;

    /// @brief Ограничивающий объем
    /// @details Не реализована по умолчанию
    virtual Box bbox() const;

    /// @brief Создать сетку общего вида
    virtual Grid make() = 0;

protected:
    /// @brief Проверить размеры сетки перед созданием
    virtual void check_size() const;

    /// @brief Проверить параметры сетки перед созданием
    virtual void check_params() const;

    /// @brief Название сеточного генератора
    std::string m_name;
};

} // namespace zephyr::mesh
