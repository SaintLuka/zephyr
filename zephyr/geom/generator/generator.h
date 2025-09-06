#pragma once

#include <string>
#include <memory>

namespace zephyr::utils { class Json; }
namespace zephyr::mesh { class AmrCells; }

namespace zephyr::geom {

struct Box;  ///< Ограничивающий объем
class Grid; ///< Сетка общего вида, создается генератором

/// @brief Базовый класс для сеточных генераторов. Содержит набор
/// виртуальных функций для создания AmrStorage с сеткой.
class Generator {
protected:
    using Json = utils::Json;
public:
    using Ptr = std::shared_ptr<Generator>;
    using Ref = const std::shared_ptr<Generator>&;

    /// @brief Базовый конструктор только с именем типа сетки
    explicit Generator(const std::string &name);

    /// @brief Виртуальный деструктор (для наследования)
    virtual ~Generator() = default;

    /// @brief Создать умный указатель по файлу конфигурации
    static Generator::Ptr create(const Json& config);

    /// @brief Проверить, что можно преобразовать к наследнику
    template <class T>
    typename std::enable_if<std::is_base_of<Generator, T>::value, bool>::type
    can_cast() { return dynamic_cast<T*>(this) != nullptr; };

    /// @brief Приведение к конкретному типу
    /// @throw segfault при невозможности приведения
    /// @code
    ///     Generator::Ptr gen;
    ///     Rectangle& rect = gen->cast<Rectangle>();
    /// @endcode
    template <class T>
    std::enable_if_t<std::is_base_of_v<Generator, T>, T&>
    cast() { return *(dynamic_cast<T*>(this)); };

    /// @brief Тип сеточного генератора
    const std::string &name() const;

    /// @brief Количество ячеек сетки
    virtual int size() const = 0;

    /// @brief Ограничивающий объем
    /// @details Не реализована по умолчанию
    virtual Box bbox() const;

    /// @brief Создать сетку общего вида
    virtual Grid make() = 0;

    /// @brief Инициализация SoA-хранилища сетки
    virtual void initialize(mesh::AmrCells& cells) { throw std::runtime_error("not implemented"); }

protected:
    /// @brief Проверить размеры сетки перед созданием
    virtual void check_size() const;

    /// @brief Проверить параметры сетки перед созданием
    virtual void check_params() const;

    /// @brief Название сеточного генератора
    std::string m_name;
};

} // namespace zephyr::mesh
