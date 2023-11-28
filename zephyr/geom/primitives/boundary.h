#pragma once

#include <string>
#include <iostream>
#include <algorithm>

namespace zephyr::geom {

/// @brief Тип граничных условий
enum class Boundary : int {
    UNDEFINED = 0,  ///< Не определено.
    ORDINARY  = 1,  ///< Обычное соседство (внутри сетки).
    WALL      = 2,  ///< Непроницаемая стенка.
    ZOE       = 3,  ///< Отображение данных запрашиваемой ячейки (простой снос).
    PERIODIC  = 4,  ///< Периодичность.
    FUNCTION  = 5,  ///< Кастомная
};

/// @brief Преобразовать тип граничного условия в строку
inline std::string to_string(Boundary type) {
    switch (type) {
        case Boundary::ORDINARY:
            return "ordinary";
        case Boundary::WALL:
            return "wall";
        case Boundary::ZOE:
            return "zoe";
        case Boundary::PERIODIC:
            return "periodic";
        case Boundary::FUNCTION:
            return "function";
        default:
            return "undefined";
    }
}

/// @brief Вывод типа гран условия в поток
inline std::ostream &operator<<(std::ostream &os, Boundary boundary) {
    os << to_string(boundary);
    return os;
}

/// @brief Граничное условие по строковому названию
inline Boundary boundary_from_string(std::string flag) {
    std::transform(flag.begin(), flag.end(), flag.begin(), ::tolower);

    if (flag == "ordinary") {
        return Boundary::ORDINARY;
    } else if (flag == "wall") {
        return Boundary::WALL;
    } else if (flag == "zoe") {
        return Boundary::ZOE;
    } else if (flag == "periodic") {
        return Boundary::PERIODIC;
    } else if (flag == "function") {
        return Boundary::FUNCTION;
    } else if ("undefined") {
        return Boundary::UNDEFINED;
    } else {
        throw std::runtime_error("Unknown boundary flag '" + flag + "'");
    }
}

} // namespace zephyr::geom