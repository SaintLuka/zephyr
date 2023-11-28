/// @file Самые базовые функции для работы с адаптивными ячейками

#pragma once

namespace zephyr::geom {

/// @brief Число вершин простой грани.
inline constexpr int VpF(int dim) { return dim < 3 ? 2 : 4; }

/// @brief Число граней у простой ячейки.
inline constexpr int FpC(int dim) { return dim < 3 ? 4 : 6; }

/// @brief Количество дочерних ячеек.
inline constexpr int CpC(int dim) { return dim < 3 ? 4 : 8; }

/// @brief Количество подграней.
inline constexpr int FpF(int dim) { return dim < 3 ? 2 : 4; }

} // namespace zephyr::geom