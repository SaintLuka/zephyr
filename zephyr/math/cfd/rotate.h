#pragma once

#include <zephyr/geom/vector.h>

namespace zephyr::math {

/// @struct Rotate
/// @brief Класс со статическими методами поворотов.
/// @details
/// Набор функций to_local поворачивает систему координат так, чтобы нормаль
/// совпала с осью Ox.
/// Набор функций to_global поворачивает систему координат так, чтобы ось Ox
/// совпала с нормалью.
/// Функции to_global и to_local являются взаимообратными.
/// Более абстрактно: функции to_local позволяют преобразовать вектор к локальной
/// системе координат, связанной с гранью ячейки с нормалью normal.
struct Rotate {
public:
    /// @brief Повернуть вектор скорости
    static void to_local(geom::Vector3d& velocity, const geom::Vector3d& normal);

    /// @brief Повернуть вектор скорости
    static void to_global(geom::Vector3d& velocity, const geom::Vector3d& normal);

    /// @brief Повернуть тензор напряжений
    static void to_local(geom::Matrix3d& stress, const geom::Vector3d& normal);

    /// @brief Повернуть тензор напряжений
    static void to_global(geom::Matrix3d& stress, const geom::Vector3d& normal);
};

}