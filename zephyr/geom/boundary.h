#pragma once

#include <string>
#include <algorithm>

namespace zephyr::geom {

/// @brief Тип граничных условий
enum class Boundary : int {
    UNDEFINED =-1,  ///< Не определено.
    INNER     = 0,  ///< Обычное соседство (внутри сетки).
    WALL      = 1,  ///< Непроницаемая стенка.
    ZOE       = 2,  ///< Отображение данных запрашиваемой ячейки (простой снос).
    PERIODIC  = 3,  ///< Периодичность.
    FUNCTION  = 4,  ///< Задаётся пользователем.
};

/// @brief Преобразовать тип граничного условия в строку
constexpr std::string to_string(Boundary type) {
    switch (type) {
        case Boundary::INNER:    return "inner";
        case Boundary::WALL:     return "wall";
        case Boundary::ZOE:      return "zoe";
        case Boundary::PERIODIC: return "periodic";
        case Boundary::FUNCTION: return "function";
        default: return "undefined";
    }
}

/// @brief Вывод типа граничного условия в консоль
inline std::ostream &operator<<(std::ostream &os, Boundary boundary) {
    os << to_string(boundary);
    return os;
}

/// @brief Граничное условие по строковому названию
inline Boundary boundary_from_string(std::string flag) {
    std::ranges::transform(flag, flag.begin(), ::tolower);

    if (flag == "inner") {
        return Boundary::INNER;
    } else if (flag == "wall") {
        return Boundary::WALL;
    } else if (flag == "zoe") {
        return Boundary::ZOE;
    } else if (flag == "periodic") {
        return Boundary::PERIODIC;
    } else if (flag == "function") {
        return Boundary::FUNCTION;
    } else {
        return Boundary::UNDEFINED;
    }
}

} // namespace zephyr::geom