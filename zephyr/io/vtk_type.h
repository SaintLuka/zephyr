#pragma once

#include <string>
#include <iostream>

namespace zephyr::io {

/// @brief Тип скаляра при сохранении в VTU, работает как enum class.
class VtkType {
public:
    enum ValueType : int {
        Undefined = 0,
        Int8 = 1, Int16 = 2, Int32 = 3, Int64 = 4,
        UInt8 = 5, UInt16 = 6, UInt32 = 7, UInt64 = 8,
        Float32 = 9, Float64 = 10
    };

    /// @brief По умолчанию Undefined
    VtkType();

    /// @brief Позволяет делать присваивание вида VtkType type = VtkType::Int32
    VtkType(const ValueType& value);

    /// @brief Позволяет делать присваивание вида VtkType type = "Int32"
    VtkType(const char* value);

    /// @brief Позволяет делать присваивание вида VtkType type = "Int32"
    VtkType(const std::string& value);

    /// @brief Позволяет делать присваивание вида VtkType type = "Int32"
    operator ValueType() const {
        return value_type_conversion();
    }

    /// @brief Конвертация в сторку C типа
    operator const char*() const {
        return const_char_conversion();
    }

    /// @brief Конвертация в C++ строку
    operator std::string() const {
        return string_conversion();
    }

    /// @brief Оператор сравнения
    bool operator==(const VtkType& value) const;

    /// @brief Оператор сравнения
    bool operator!=(const VtkType& value) const;

    /// @brief Определен ли тип
    bool is_undefined() const;

    /// @brief Является ли тип целым числом
    bool is_signed() const;

    /// @brief Явялется ли тип натуральным числом
    bool is_unsigned() const;

    /// @brief Является ли тип значением с плавающей точкой
    bool is_floating() const;

    /// @brief Размер типа в байтах (sizeof типа)
    size_t size() const;

    /// @brief Конвертация в VtkType из типа C++.
    /// Пример. VtkType type = VtkType::get<int>();
    template <class T>
    static VtkType get() {
        if (!std::is_arithmetic<T>::value) {
            return VtkType::Undefined;
        }
        if (std::is_floating_point<T>::value) {
            return get_floating(sizeof(T));
        }
        if (std::is_signed<T>::value) {
            return get_signed(sizeof(T));
        }
        return get_unsigned(sizeof(T));
    }

private:
    const char* const_char_conversion() const;
    ValueType value_type_conversion() const;
    std::string string_conversion() const;

    static VtkType get_signed(size_t size);
    static VtkType get_unsigned(size_t size);
    static VtkType get_floating(size_t size);

    ValueType m_value;
};

/// @brief Оператор для вывода типа в консоль
std::ostream& operator<<(std::ostream& os, const VtkType& vtk_type);

} // namespace zephyr::io

