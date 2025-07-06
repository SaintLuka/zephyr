#pragma once

#include <string>
#include <iostream>
#include <variant>
#include <functional>

#include <zephyr/geom/vector.h>
#include <zephyr/io/vtk_type.h>

// Forward declaration
namespace zephyr::mesh  { class EuCell; }
namespace zephyr::utils { template<typename T> struct Storable; }


namespace zephyr::io {

/// @brief Тип функции для записи переменных, позволяет инициализировать
/// функцию записи переменных через лямбда функцию. Далее пример использования.
/// позволяет сократить размер выходного файла.
/// @code
/// WriteFunction<float> ev_energy =
///     [](EuCell& cell, float* out) {
///         out[0] = static_cast<float>(cell.energy / 1.6e-19);
///     };
/// @endcode
template <typename T>
using WriteSoaCell = std::function<void(mesh::EuCell&, T*)>;

/// @brief Класс для записи переменных в VTU файл, каждой переменной для
/// записи должен соответствовать экземпляр Variable.
class Variable {
public:

    /// @brief Создание дескриптора без описания запрещено
    Variable() = delete;

    /// @brief Создание дескриптора по имени.
    /// @name Имя переменной
    /// @details Функция актуальна для некоторых предопеределенных имен:
    /// "coords", "center", "volume"...
    explicit Variable(const char *name);

    /// @brief Аналогично конструктору Variable(const char* )
    explicit Variable(const std::string &name);

    /// @brief Создать переменную с полным описанием
    /// @param name Имя переменной
    /// @param n_components Размер вектора для хранения переменной
    /// @param func Функция записи переменной по заданному указателю
    /// @tparam T Тип шаблона важен для правильной дедукции VtkType
    /// Следующий код позволяет добавить векторную переменную "momentum",
    /// которая позволяет записывать две компоненты импульса в формате
    /// Float32. В целом, использование чисел с одинарной точностью
    /// позволяет сократить размер выходного файла.
    /// @code
    /// Variable fd("momentum", 2,
    ///     WriteFunction<float>([](AmrStorage::Item cell, float* out) {
    ///         out[0] = static_cast<float>(cell.mass * cell.velocity.x);
    ///         out[1] = static_cast<float>(cell.mass * cell.velocity.y);
    ///     }));
    /// @endcode
    template<class T>
    Variable(const char *name,
             int n_components,
             const WriteSoaCell<T> &func) {

        m_name = name;
        m_type = VtkType::get<T>();
        m_n_components = n_components;
        m_write = [func](mesh::EuCell &cell, void *out) {
            func(cell, (T *) out);
        };
    }

    /// Тип vtk_type выводится из T.
    template<class T>
    Variable(const std::string &name, int n_components, const WriteSoaCell<T> &func)
            : Variable(name.c_str(), n_components, func) {}

    /// @brief Имя переменной
    std::string name() const;

    /// @brief Переменная эйлеровой ячейки
    bool is_eu_cell() const;

    /// @brief Тип переменной
    VtkType type() const;

    /// @brief Число компонент для векторной переменной
    int n_components() const;

    /// @brief Является ли переменная скаляром
    bool is_scalar() const;

    /// @brief Размер переменной в байтах (аналог sizeof)
    size_t size() const;

    /// @brief Основная функция класса. Запись переменной из ячейки в поток.
    void write(mesh::EuCell &cell, void *out) const;

private:
    std::string m_name;
    VtkType m_type;
    int m_n_components;

    // Три функции записи для разных типов элементов,
    // в каждом экземпляре актуальна только одна.
    WriteSoaCell<void> m_write = nullptr;
};

} // namespace zephyr::io

