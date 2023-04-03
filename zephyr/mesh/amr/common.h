/// @file Файл содержит набор определений и реализацию некоторых простых функций,
/// используемых в алгоритмах адаптации.
/// Данный файл не устанавливается при установке zephyr, все изложенные описания
/// алгоритмов и комментарии к функциям предназначены исключительно для разработчиков.

#pragma once

#define SCRUTINY 1          ///< Тщательная проверка вычислений
#define CHECK_PERFORMANCE 1 ///< Выводить производительность частей алгоритма
#define FAST_BALANCING 1    ///< Быстрая функция балансировки флагов

#if SCRUTINY
/// @brief Бросает исключение, если условие не выполняется
#define scrutiny_check(condition, message) if (!(condition)) throw std::runtime_error(message);
#else
/// @brief Опция SCRUTINY отключена, ничего не происходит
#define scrutiny_check(...)
#endif

#include <array>
#include <zephyr/geom/base.h>

namespace zephyr { namespace mesh { namespace amr {

using namespace zephyr::geom;

/// @brief Локальные индексы ячеек, прилегающих к стороне.
template<int dim>
inline constexpr std::array<std::array<int, VpF(dim)>, FpC(dim)> get_children_by_side();

/// @brief Локальные индексы ячеек, прилегающих к стороне (2D).
template<>
inline constexpr std::array<std::array<int, 2>, 4> get_children_by_side<2>() {
    return {
            std::array<int, 2>({0, 2}),
            std::array<int, 2>({1, 3}),
            std::array<int, 2>({0, 1}),
            std::array<int, 2>({2, 3})
    };
}

/// @brief Локальные индексы ячеек, прилегающих к стороне (3D).
template<>
inline constexpr std::array<std::array<int, 4>, 6> get_children_by_side<3>() {
    return {
            std::array<int, 4>({0, 2, 4, 6}),
            std::array<int, 4>({1, 3, 5, 7}),
            std::array<int, 4>({0, 1, 4, 5}),
            std::array<int, 4>({2, 3, 6, 7}),
            std::array<int, 4>({0, 1, 2, 3}),
            std::array<int, 4>({4, 5, 6, 7})
    };
}

/// @brief Индексы подграней на стороне в списке faces.
template<int dim>
inline std::array<int, FpF(dim)> subface_sides(int side);

/// @brief  Индексы подграней на стороне в списке faces (2D).
template<>
inline std::array<int, 2> subface_sides<2>(int side) {
    return {side, side + 6};
}

/// @brief Индексы подграней на стороне в списке faces (3D).
template<>
inline std::array<int, 4> subface_sides<3>(int side) {
    return {side, side + 6, side + 12, side + 18};
}

} // namespace amr
} // namespace mesh
} // namespace zephyr
