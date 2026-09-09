#pragma once

#include <string>
#include <functional>
#include <variant>

#include <zephyr/io/vtk_type.h>

// Forward declaration
namespace zephyr::mesh {
template<typename T>
class Storable;
class EuCell;
class EuNode;
}

namespace zephyr::io {

/// @brief Тип функции для записи переменных, позволяет инициализировать
/// функцию записи переменных через лямбда функцию. Далее пример использования.
/// Позволяет сократить размер выходного файла.
/// @code
/// WriteFunction<float> ev_energy =
///     [](EuCell& cell, float* out) {
///         out[0] = static_cast<float>(cell.energy / 1.6e-19);
///     };
/// @endcode
template <typename T>
using WriteCell = std::function<void(mesh::EuCell&, T*)>;

/// @brief Тип функции для записи переменных из узлов
template <typename T>
using WriteNode = std::function<void(mesh::EuNode&, T*)>;

/// @brief Класс для записи переменных в VTU файл, каждой переменной для
/// записи должен соответствовать экземпляр Variable.
class Variable {
public:
    /// @brief Создание дескриптора без описания запрещено
    Variable() = delete;

    /// @brief Создание дескриптора по имени.
    /// @name Имя переменной
    /// @details Функция актуальна для некоторых предопределенных имен:
    /// "coord", "center", "volume"...
    explicit Variable(std::string_view name);

    /// @brief Создать переменную с полным описанием
    /// @param name Имя переменной
    /// @param n_components Размер вектора для хранения переменной
    /// @param func Функция записи переменной по заданному указателю
    /// @tparam T Тип шаблона важен для правильной дедукции VtkType
    ///
    /// Следующий код добавляет векторную переменную "momentum", которая
    /// позволяет записывать две компоненты импульса в формате Float32.
    /// @code
    /// Variable fd("momentum", 2,
    ///     WriteFunction<float>([](EuCell& cell, float* out) {
    ///         out[0] = static_cast<float>(cell.mass * cell.velocity.x);
    ///         out[1] = static_cast<float>(cell.mass * cell.velocity.y);
    ///     }));
    /// @endcode
    template<class T>
    Variable(std::string_view name, int n_components, const WriteCell<T> &func) {
        name_ = name;
        type_ = VtkType::get<T>();
        n_components_ = n_components;
        write_ = [func](mesh::EuCell &cell, void *out) {
            func(cell, static_cast<T *>(out));
        };
    }

    template<class T>
    Variable(std::string_view name, int n_components, const WriteNode<T> &func) {
        name_ = name;
        type_ = VtkType::get<T>();
        n_components_ = n_components;
        write_ = [func](mesh::EuNode &node, void *out) {
            func(node, static_cast<T *>(out));
        };
    }

    /// @brief Имя переменной
    std::string name() const { return name_; }

    /// @brief Тип переменной
    VtkType type() const { return type_; }

    /// @brief Число компонент для векторной переменной
    int n_components() const { return n_components_; }

    /// @brief Является ли переменная скаляром
    bool is_scalar() const { return n_components_ < 2; }

    /// @brief Размер переменной в байтах (аналог sizeof)
    size_t size() const { return n_components_ * type_.size(); }

    /// @brief Переменная для записи сеточных данных?
    bool cell_data() const;

    /// @brief Переменная для записи сеточных данных?
    bool node_data() const;

    /// @brief Основная функция класса. Запись переменной из ячейки в буфер.
    void write(mesh::EuCell &cell, void *out) const;

    /// @brief Основная функция класса. Запись переменной из узла в буфер.
    void write(mesh::EuNode &node, void *out) const;

private:
    std::string name_;  ///< Имя переменной
    VtkType type_;      ///< Тип переменной
    int n_components_;  ///< Число компонент (для вектора)

    /// @brief Функция записи
    std::variant<WriteCell<void>, WriteNode<void>> write_ = {};
};

} // namespace zephyr::io

