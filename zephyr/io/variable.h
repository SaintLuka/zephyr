#pragma once

#include <string>
#include <iostream>
#include <functional>

#include <zephyr/io/vtk_type.h>

// Forward declaration
namespace zephyr::mesh {
template<typename T>
struct Storable;
class EuCell;
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

/// @brief Класс для записи переменных в VTU файл, каждой переменной для
/// записи должен соответствовать экземпляр Variable.
class Variable {
public:
    /// @brief Создание дескриптора без описания запрещено
    Variable() = delete;

    /// @brief Создание дескриптора по имени.
    /// @name Имя переменной
    /// @details Функция актуальна для некоторых предопределенных имен:
    /// "coords", "center", "volume"...
    explicit Variable(const char *name);

    /// @brief Аналогично конструктору Variable(const char* )
    explicit Variable(const std::string &name);

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
    Variable(const char *name, int n_components, const WriteCell<T> &func) {
        m_name = name;
        m_type = VtkType::get<T>();
        m_n_components = n_components;
        m_write = [func](mesh::EuCell &cell, void *out) {
            func(cell, static_cast<T *>(out));
        };
    }

    /// Тип VtkType выводится из T.
    template<class T>
    Variable(const std::string &name, int n_components, const WriteCell<T> &func)
            : Variable(name.c_str(), n_components, func) { }

    /// @brief Имя переменной
    std::string name() const { return m_name; }

    /// @brief Тип переменной
    VtkType type() const { return m_type; }

    /// @brief Число компонент для векторной переменной
    int n_components() const { return m_n_components; }

    /// @brief Является ли переменная скаляром
    bool is_scalar() const { return m_n_components < 2; }

    /// @brief Размер переменной в байтах (аналог sizeof)
    size_t size() const { return m_n_components * m_type.size(); }

    /// @brief Основная функция класса. Запись переменной из ячейки в поток.
    void write(mesh::EuCell &cell, void *out) const;

private:
    std::string m_name;  ///< Имя переменной
    VtkType m_type;      ///< Тип переменной
    int m_n_components;  ///< Число компонент (для вектора)

    WriteCell<void> m_write = nullptr; ///< Функция записи
};

} // namespace zephyr::io

