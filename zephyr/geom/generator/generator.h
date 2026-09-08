#pragma once

#include <string>
#include <memory>
#include <zephyr/mesh/euler/amr_cells.h>

namespace zephyr::utils { class Json; }
namespace zephyr::mesh { class AmrCells; }

namespace zephyr::geom {

struct Box;  ///< Ограничивающий объем
class Grid;  ///< Сетка общего вида, создается генератором

/// @brief Максимальное число ячеек сетки
constexpr int max_grid_size = 200'000'000;

/// @brief Базовый класс для сеточных генераторов.
/// Виртуальная функция make() -> Grid создает сетку.
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

    /// @brief Создать умный указатель по конфигу
    static Generator::Ptr create(const Json& config);

    /// @brief Проверить, что можно преобразовать к наследнику
    template <class T>
    std::enable_if_t<std::is_base_of_v<Generator, T>, bool>
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

    /// @brief Структурированная сетка (в редких случаях true)
    virtual bool structured() const { return false; }

    /// @brief Сетка с осевой симметрией?
    bool axial() const { return axial_; }

    /// @brief Адаптивная сетка?
    bool adaptive() const { return adaptive_; }

    /// @brief Линейная адаптивная сетка?
    bool linear() const { return linear_; }

    /// @brief Использовать осевую симметрию
    virtual void set_axial(bool axial);

    /// @brief Использовать адаптацию
    virtual void set_adaptive(bool adaptive);

    /// @brief Использовать линейную адаптацию
    virtual void set_linear(bool linear);

    /// @brief Ограничивающий объем
    /// @details Не реализована по умолчанию
    virtual Box bbox() const;

    /// @brief Создать сетку общего вида
    virtual Grid make() const = 0;

    /// @brief Некоторые генераторы могут напрямую инициализировать хранилище
    virtual bool can_make_cells() const { return false; }

    /// @brief Инициализация SoA-хранилища сетки
    virtual mesh::AmrCells make_cells(bool unique_nodes) const {
        throw std::runtime_error("Generator::make_cells: not implemented");
    }

protected:
    /// @brief Проверить размеры сетки перед созданием
    virtual void check_size(size_t size) const;

    std::string name_;     ///< Название сеточного генератора
    bool axial_{false};    ///< Сетка с осевой симметрией
    bool adaptive_{false}; ///< Адаптивная сета
    bool linear_{true};    ///< Простая адаптивная сетка
};

} // namespace zephyr::mesh
