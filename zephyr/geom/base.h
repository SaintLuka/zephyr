/// @file Самые базовые функции для работы с адаптивными ячейками

#pragma once

namespace zephyr { namespace geom {

/// @brief Число вершин простой грани.
inline constexpr int VpF(int dim) { return dim < 3 ? 2 : 4; }

/// @brief Число вершин простой ячейки.
inline constexpr int VpC(int dim) { return dim < 3 ? 4 : 8; }

/// @brief Число граней у простой ячейки.
inline constexpr int FpC(int dim) { return dim < 3 ? 4 : 6; }

/// @brief Количество дочерних ячеек.
inline constexpr int CpC(int dim) { return dim < 3 ? 4 : 8; }

/// @brief Количество подграней.
inline constexpr int FpF(int dim) { return dim < 3 ? 2 : 4; }

/*
    Далее несколько функций, упрощающих индексацию
    в четырехугольных ячейках и шестигранниках, особенно при использовании AMR.
    Отображения строятся из множества трехмерных индексов на одномерный индекс.
    Коротким (s/short) индексом будем называть трехмерный индекс из
    множества {0, 1}^3 или одномерный индекс из [0, 8).
    Длинным (w/wide) индексом будем называть трехмерный индекс из
    множества {0, 1, 2}^3 или одномерный индекс из [0, 27).
 */

/// @brief Отображение {0, 1}^3 -> [0, 8)
/// @return Индекс [0, 8) вершины ячейки.
inline constexpr int iss(int i, int j, int k = 0) {
    return 4 * k + 2 * j + i;
}

/// @brief Отображение {0, 2}^3 -> [0, 8)
/// @return Индекс [0, 8) вершины ячейки.
inline constexpr int iws(int i, int j, int k = 0) {
    return 2 * k + j + i / 2;
}

/// @brief Отображение {0, 1}^3 -> [0, 27)
/// @return Индекс [0, 27) вершины ячейки.
inline constexpr int isw(int i, int j, int k = 0) {
    return 2 * (9 * k + 3 * j + i);
}

/// @brief Отображение {0, 1, 2}^3 -> [0, 27)
/// @return Индекс [0, 27) вершины ячейки.
inline constexpr int iww(int i, int j, int k = 0) {
    return 9 * k + 3 * j + i;
}

/// @brief Динамическое отображение. Отображает короткий индекс {0, 1}^3
/// на одномерный индекс, соответствующий type::_vertices_::max_size
inline constexpr int iv(int i, int j, int k = 0) {
    return isw(i, j, k);
}

} // geom
} // zephyr